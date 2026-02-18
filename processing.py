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
        self.repeat_motif: Dict[Tuple[str, int, int], List[str]] = {}
        self.walks: Dict[Tuple[str, str], Dict[Tuple[int, int], List[int]]] = {}
        self.repeat_nodes: Dict[Tuple[str, int, int], Set[int]] = {}
        self.flanking_nodes: Dict[Tuple[str, int, int], Tuple[int, int, int, int]] = {}
        self.pangenome_reading(graph_file)

        def complementary(char: str) -> str:
            match char:
                case 'G':
                    return 'C'
                case 'C':
                    return 'G'
                case 'A':
                    return 'T'
                case 'T':
                    return 'A'
                case _:
                    return ''

        canonical = "TTAGGG"
        self.telomeres: Set[str] = {canonical[i:] + canonical[:i] for i in range(6)}
        rev_canonical: str = "".join(map(complementary, reversed(canonical)))
        self.telomeres.update(rev_canonical[i:] + rev_canonical[:i] for i in range(6))
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
    
    def binary_search_nodes(self, start: int, end: int) -> Tuple[Set[int], int, int]:
        len_l = len(self.nodes_ref_list)
        l, r = 0, len_l - 1
        output_nodes: Set[int] = set()
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
        i_iter_r: int = i + 1
        output_nodes.add(self.nodes_ref_list[i][0])
        while i_iter_r < len_l and self._overlaps(start, end, i_iter_r): # then I iterate over its neighbors to construct a sequence of nodes covering it
            output_nodes.add(self.nodes_ref_list[i_iter_r][0])
            i_iter_r += 1
        i_iter_l: int = i - 1
        while i_iter_l > -1 and self._overlaps(start, end, i_iter_l):
            output_nodes.add(self.nodes_ref_list[i_iter_l][0])
            i_iter_l -= 1
        return output_nodes, self.nodes_ref_list[i_iter_l+1][0], self.nodes_ref_list[i_iter_r-1][0]


    def tandem_reading(self, repeat_file: str) -> None:

        def is_telomere_motif(motif: str) -> bool:
            if len(motif) % 6 != 0:
                return False
            elif len(motif) == 12:
                return motif[:6] in self.telomeres and motif[6:] in self.telomeres
            else:
                return motif in self.telomeres

        self.r = repeat_file
        file = open(repeat_file, 'r')
        mem: List[Tuple[str, int, int, List[str]]] = []
        for line in file.readlines():
            contig, start, end, num, _, seq, _ = line.split()
            if len(seq) > 12 or is_telomere_motif(seq):
                continue
            start, end = map(int, (start, end))
            if mem and end - mem[-1][2] <= 50 and contig == mem[-1][0]:
                p_contig, p_start, p_end, p_seq_l = mem.pop()
                mem.append((p_contig, p_start, end, p_seq_l + [seq]))
                continue
            mem.append((contig, start, end, [seq]))
        for contig, start, end, seq_l in mem:
            nodes_sequence, left_flanking, right_flanking = self.binary_search_nodes(start, end)
            if len(nodes_sequence) > 1:
                self.repeat_nodes[contig, start, end] = nodes_sequence
                self.repeat_motif[contig, start, end] = seq_l
                left_split = start - self.nodes_ref_ind[left_flanking]
                right_split = end - self.nodes_ref_ind[right_flanking]
                self.flanking_nodes[contig, start, end] = (left_flanking, left_split, right_flanking, right_split)
        mem.clear()
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
            paths_seqs: Dict[Tuple[str, int, int], Dict[str, Tuple[str, str, str]]],
            contig: str,
            divisions: Dict[Tuple[str, int, int], Dict[str, Set[int]]]) -> None:
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
        node_spans: List[Tuple[int, int, int]] = []
        cursor: int = 0
        seq: str = ""
        for x in stack:
            s = rev_if_neg(x)
            seq += s
            node_spans.append((abs(x), cursor, cursor + len(s)))
            cursor += len(s)
        current_start: int = current
        current_end: int = current_start + len(seq)
        flanks: Tuple[int, int, int, int] = self.flanking_nodes[ref_key]
        left, repeat, right = self.split_by_flanks(seq, node_spans, flanks)
        len_left = len(left)
        local_divisions: Set[int] = {x[1] - len_left for x in node_spans if len_left < x[1]}
        #if len(local_divisions) <= 1:
            #print(" ".join(map(lambda x : self.nodes[abs(x)], stack)))
        i_beg = i - 1
        while len_left < 50 and i_beg >= 0:
            len_left += len(self.nodes[abs(walk[i_beg])])
            i_beg  -= 1
        if i_beg < i:
            left = ("".join(map(rev_if_neg, walk[i_beg:i]))+left)[-50:]
        right_len = len(right)
        i_end = i + len(stack) + 1
        while right_len < 50 and i_end < len(walk):
            right_len += len(self.nodes[abs(walk[i_end])])
            i_end += 1
        if i_end > i + len(stack):
            right = (right+"".join(map(rev_if_neg, walk[i+len(stack):i_end])))[:50]

        if (ref_contig, start, end) not in paths:
            paths[ref_contig, start, end] = {contig: (current_start, current_end)}
            paths_seqs[ref_contig, start, end] = {contig: (left, repeat, right)}
            divisions[ref_contig, start, end] = {contig: local_divisions}
        elif contig not in paths[ref_contig, start, end] or (paths[ref_contig, start, end][contig][1] -
                paths[ref_contig, start, end][contig][0]) < current_end - current_start:
            paths[ref_contig, start, end][contig] = (current_start, current_end)
            paths_seqs[ref_contig, start, end][contig] = (left, repeat, right)
            divisions[ref_contig, start, end][contig] = local_divisions


    def analysis(
            self,
            mismatch_penalty: int,
            match_score: int,
            threshold: int
    ):
        self.nodes_ref_ind.clear()
        self.nodes_ref_list.clear()
        node_occurrences = defaultdict(list)
        for (_, contig), nodes_dict in self.walks.items():
            for (start_w, _), walk in nodes_dict.items():
                for i, node in enumerate(walk):
                    node_occurrences[node].append((walk, i, start_w, contig))
        paths: Dict[Tuple[str, int, int], Dict[str, Tuple[int, int]]] = {}
        paths_seqs: Dict[Tuple[str, int, int], Dict[str, Tuple[str, str, str]]] = {}
        divisions: Dict[Tuple[str, int, int], Dict[str, Set[int]]] = {}
        for ref_key, nodes_rf in self.repeat_nodes.items():
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
                        divisions=divisions
                    )
        with open(self.l, 'w') as log:
            for key, info in paths_seqs.items():
                print(
                    f"Tandem repeat motif 5'-{"-*-".join(self.repeat_motif[key])}"
                    f"-3'\nIn {key[0]} at [ {key[1]}, {key[2]} ]:",
                    file=log
                )
                blocks: List[Tuple[str, str, str, str]] = []
                ref_ind: Optional[int] = None
                ind: int = 0
                for contig, (left, repeat, right) in info.items():
                    if contig == key[0]:
                        ref_ind = ind
                    blocks.append((left, repeat, right, contig))
                    ind += 1
                del ind
                blocks[0], blocks[ref_ind] = blocks[ref_ind], blocks[0]
                self.pretty_print(blocks, key, paths, divisions, log=log)

    def split_by_flanks(self, seq: str, node_spans, flanks) -> Tuple[str, str, str]:
        left_node, left_off, right_node, right_off = flanks
        repeat_start = None
        repeat_end = None

        for node, s, e in node_spans:
            if node == left_node and repeat_start is None:
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
        return left[:50] , repeat, right[:50]

    def pretty_print(self, blocks: List[Tuple[str, str, str, str]], ref_key: Tuple[str, int, int],
                     paths: Dict[Tuple[str, int, int], Dict[str, Tuple[int, int]]],
                     divisions: Dict[Tuple[str, int, int], Dict[str, Set[int]]], log):
        max_rep: int = 0
        for n, (l, t, r, contig) in enumerate(blocks):
            new_t = ""
            for m, x in enumerate(t):
                if m in divisions[ref_key][contig]:
                    new_t += '_'
                new_t += x
            if len(new_t) > max_rep:
                max_rep = len(new_t)
            blocks[n] = (l, new_t, r, contig)
        max_left: int = max(len(l) for l, _, _, _ in blocks)
        ref = ref_key[0]
        for left, rep, right, contig in blocks:
            approx_start: int = paths[ref_key][contig][0]
            approx_end: int = paths[ref_key][contig][1]
            if contig != ref:
                print(f"Found in {contig} at [ {approx_start}, {approx_end} ]:", file=log)
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
