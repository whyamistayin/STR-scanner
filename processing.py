import argparse
from typing import Tuple, List, Dict, Set, Optional
from collections import defaultdict

class Processing:

    def __init__(self, graph_file: str, repeat_file: str, logfile: str) -> None:
        self.l = logfile
        self.g = graph_file
        self.r = repeat_file
        self.nodes_ref_list: List[Tuple[int, str, int]] = []
        self.nodes_ref_ind: Dict[int, int] = {}
        self.nodes: Dict[int, str] = {}
        self.walks: Dict[Tuple[str, str], Dict[Tuple[int, int], List[int]]] = {}
        self.repeat_nodes: Dict[Tuple[str, int, int], Set[int]] = {}
        self.flanking_nodes: Dict[Tuple[str, int, int], Tuple[int, int, int, int]] = {}
        self.pangenome_reading(graph_file)
        self.tandem_reading(repeat_file)

    def __str__(self) -> str:
        return " | ".join([self.g, self.r, self.l])
    
    __repr__ = __str__

    def _parse_walk(self, walk: str) -> List[int]:
        state: int = 0 # 0 is determining direction, 1 reads number
        direction: int = 1
        read: str = ""
        output: List[int] = []
        for c in walk:
            if state == 1:
                if c in ['>', '<']:
                    state = 0
                else:
                    read += c
            if state == 0:
                if read:
                    output.append(int(read) * direction)
                    read = ""
                direction = -1 if c == '<' else 1
                state = 1
        return output

    def pangenome_reading(self, graph_file: str) -> None:
        graph = open(graph_file, 'r')
        for line in graph.readlines():
            if line.startswith("S"):
                line = line.split()
                id, contains = line[1:3]
                id = int(id)
                if len(line) > 3:
                    index: int = int(line[4].split(':')[-1])
                    self.nodes_ref_list.append((id, contains, index))
                    self.nodes_ref_ind[id] = index
                self.nodes[id] = contains
            elif line.startswith("W"):
                ref, _, con, st, end, walk = line.split()[1:]
                processed = self._parse_walk(walk)
                #print(processed)
                if (ref, con) not in self.walks:
                    self.walks[ref, con] = {(int(st), int(end)): processed[:]}
                else:
                    self.walks[ref, con][int(st), int(end)] = processed[:]
        graph.close()

    def _overlaps(self, start: int, end: int, i: int) -> bool:
        if not (0 <= i < len(self.nodes_ref_list)):
            return False
        st_node = self.nodes_ref_list[i][2]
        end_node = st_node + len(self.nodes_ref_list[i][1]) - 1
        return start <= end_node and st_node <= end
    
    def binary_search_nodes(self, start: int, end: int) -> Tuple[Set[int], Tuple[int, int]]:
        len_l = len(self.nodes_ref_list)
        l, r = 0, len_l - 1
        output_nodes: List[int] = []
        while l <= r: # I am searching for at least some node inside a tandem repeat
            i = (r + l) // 2
            if self._overlaps(start, end, i):
                break
            if self.nodes_ref_list[i][2] < start:
                l = i + 1
            elif self.nodes_ref_list[i][2] > end:
                r = i - 1
            else:
                raise IndexError()
        i_iter: int = i + 1
        output_nodes.append(self.nodes_ref_list[i][0])
        while i_iter < len_l and self._overlaps(start, end, i_iter): # then I iterate over its neighbors to construct a sequence of nodes covering it
            output_nodes.append(self.nodes_ref_list[i_iter][0])
            i_iter += 1
        i_iter: int = i - 1
        output_nodes2 = []
        while i_iter > -1 and self._overlaps(start, end, i_iter):
            output_nodes2.append(self.nodes_ref_list[i_iter][0])
            i_iter -= 1
        return set(output_nodes+output_nodes2), (output_nodes2[-1] if output_nodes2 else output_nodes[0],
                                                      output_nodes[-1])


    def tandem_reading(self, repeat_file: str) -> None:
        self.r = repeat_file
        file = open(repeat_file, 'r')
        for line in file.readlines():
            contig, start, end, num, _, seq, _ = line.split()
            start, end = map(int, (start, end))
            nodes_sequence, (left_flanking, right_flanking) = self.binary_search_nodes(start, end)
            if len(nodes_sequence) > 1:
                self.repeat_nodes[contig, start, end] = nodes_sequence
                left_split = start - self.nodes_ref_ind[left_flanking]
                right_split = end - self.nodes_ref_ind[right_flanking]
                self.flanking_nodes[contig, start, end] = (left_flanking, left_split, right_flanking, right_split)
        file.close()

    def _evaluate(self,
            i: int,
            walk: List[int],
            threshold: int,
            match_score: int,
            mismatch_penalty: int,
            nodes_rf: Set[int],
            current: int,
            ref_key: Tuple[str, int, int],
            paths: Dict[Tuple[str, int, int], Dict[str, Tuple[int, int]]],
            paths_seqs: Dict[Tuple[str, int, int], Dict[str, Tuple[str, List[Tuple[int, int, int]]]]],
            contig: str) -> None:
        nodes = self.nodes
        rev_if_neg = lambda x: nodes[-x][::-1] if x < 0 else nodes[x]
        ref_contig, start, end = ref_key
        score: int = threshold + match_score * 2
        stack: List[int] = []
        for k in range(i, len(walk)):
            if abs(walk[k]) in nodes_rf:
                score += match_score
            else:
                score += mismatch_penalty
            if score < threshold:
                break
            stack.append(walk[k])
        while abs(stack[-1]) not in nodes_rf:
            stack.pop()
        seq_parts: List[str] = []
        node_spans: List[Tuple[int, int, int]] = []
        cursor: int = 0
        for x in stack:
            s = rev_if_neg(x)
            seq_parts.append(s)
            node_spans.append((abs(x), cursor, cursor + len(s)))
            cursor += len(s)
        seq = "".join(seq_parts)
        current_start: int = current
        current_end: int = current_start + len(seq)
        if (ref_contig, start, end) not in paths:
            paths[ref_contig, start, end] = {contig:  (current_start, current_end)}
            paths_seqs[ref_contig, start, end] = {contig: (seq, node_spans)}
        else:
            if contig not in paths[ref_contig, start, end]:
                paths[ref_contig, start, end][contig] = (current_start, current_end)
                paths_seqs[ref_contig, start, end][contig] = (seq, node_spans)
            elif (paths[ref_contig, start, end][contig][1] - paths[ref_contig, start, end][contig][0]) <\
                current_end - current_start:
                paths[ref_contig, start, end][contig] = (current_start, current_end)
                paths_seqs[ref_contig, start, end][contig] = (seq, node_spans)


    def analysis(
            self,
            mismatch_penalty: int,
            match_score: int,
            threshold: int
    ):
        node_occurrences = defaultdict(list)
        for (_, contig), nodes_dict in self.walks.items():
            for (start_w, _), walk in nodes_dict.items():
                for i, node in enumerate(walk):
                    node_occurrences[node].append((walk, i, start_w, contig))
        paths: Dict[Tuple[str, int, int], Dict[str, Tuple[int, int]]] = {}
        paths_seqs: Dict[Tuple[str, int, int], Dict[str, Tuple[str, List[Tuple[int, int, int]]]]] = {}
        for ref_key, nodes_rf in self.repeat_nodes.items():
            ref_contig, start, end = ref_key
            for node in nodes_rf:
                for walk, i, start_w, contig in node_occurrences.get(node, []):
                    self._evaluate(
                        i=i,
                        walk=walk,
                        threshold=threshold,
                        match_score=match_score,
                        mismatch_penalty=mismatch_penalty,
                        nodes_rf=nodes_rf,
                        current=start_w,
                        ref_key=ref_key,
                        paths=paths,
                        paths_seqs=paths_seqs,
                        contig=contig,
                    )
        with open(self.l, 'w') as log:
            for key, info in paths_seqs.items():
                print(
                    f"Tandem repeat in {key[0]} at [ {key[1]}, {key[2]} ]:",
                    file=log
                )

                flanks: Tuple[int, int, int, int] = self.flanking_nodes[key]
                blocks: List[Tuple[str, str, str, str]] = []
                ref_ind: Optional[int] = None

                ind: int = 0
                for contig, (seq, node_spans) in info.items():
                    split: Tuple[str, str, str] = self.split_by_flanks(seq, node_spans, flanks)
                    if contig == key[0]:
                        ref_ind = ind
                    blocks.append((*split, contig))
                    ind += 1

                blocks[0], blocks[ref_ind] = blocks[ref_ind], blocks[0]
                self.pretty_print(blocks, key, paths, log=log)

    def split_by_flanks(self, seq, node_spans, flanks):
        left_node, left_off, right_node, right_off = flanks
        repeat_start = None
        repeat_end = None

        for node, s, e in node_spans:
            if (node == left_node or self.nodes[node] == self.nodes[left_node]) and repeat_start is None:
                repeat_start = s + left_off
            if node == right_node:
                repeat_end = s + right_off
            # Fallbacks (important for graph weirdness)
        if repeat_start is None:
            #print(left_node, [x[0] for x in node_spans])
            repeat_start = 0
        if repeat_end is None:
            repeat_end = len(seq)

        left = seq[:repeat_start]
        repeat = seq[repeat_start:repeat_end]
        right = seq[repeat_end:]
        return left, repeat, right

    def pretty_print(self, blocks: List[Tuple[str, str, str, str]], ref_key: Tuple[str, int, int],
                     paths: Dict[Tuple[str, int, int], Dict[str, Tuple[int, int]]], log):
        max_left: int = max(len(l) for l, _, _, _ in blocks)
        max_rep: int = max(len(r) for _, r, _, _ in blocks)
        ref = ref_key[0]
        for left, rep, right, contig in blocks:
            true_start: int = paths[ref_key][contig][0] + len(left)
            true_end: int = paths[ref_key][contig][1] - len(right)
            if contig != ref:
                print(f"Found in {contig} at [ {true_start}, {true_end} ]:", file=log)
            print(f"{left:>{max_left}}  {rep:<{max_rep}}  {right}", file=log)
        print(file=log)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', "--graph", required=True, help="pangenome-graph")
    parser.add_argument('-r', "--repeats", required=True, help="tandem-repeats")
    parser.add_argument('-l', "--log", required=True, help="logfile")
    args = parser.parse_args()
    g, r, l = args.graph, args.repeats, args.log
    processing = Processing(g, r, l)
    #import time
    #start = time.time()
    processing.analysis(-1, 1, 2)
    #print("Computation ended in %f seconds" % (time.time() - start))
