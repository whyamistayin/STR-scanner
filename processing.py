"""
A script for mapping tandem repeat loci predictions from the reference to the graph and path reconstruction.
"""
import argparse
import numpy as np
from typing import Tuple, List, Dict, Set, Optional
from collections import defaultdict, Counter
from sys import stdin
from collections import deque

class Processing:

    def __init__(self, repeat_file: str, logfile: str, chromosome: str, reference: str) -> None:
        self.l = logfile
        self.r = repeat_file
        self.c = f"chr{chromosome}"
        self.reference = reference
        self.nodes_ref_list: List[Tuple[int, str, int]] = []
        self.nodes_ref_ind: Dict[int, int] = {}
        self.nodes: Dict[int, str] = {}
        self.repeat_motif: Dict[Tuple[str, int, int], List[str]] = {}
        self.walks: Dict[Tuple[str, str], Dict[Tuple[int, int], List[int]]] = {}
        self.repeat_nodes: Dict[Tuple[str, int, int], Set[int]] = {}
        self.flanking_nodes: Dict[Tuple[str, int, int], Tuple[int, int, int, int]] = {}
        self.pangenome_reading()

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
        return " | ".join([self.r, self.l])
    
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

    def pangenome_reading(self) -> None:
        for line in map(lambda x: x.strip(), stdin):
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

    def _overlaps(self, start: int, end: int, i: int) -> bool:
        if not (0 <= i < len(self.nodes_ref_list)):
            return False
        st_node = self.nodes_ref_list[i][2]
        end_node = st_node + len(self.nodes_ref_list[i][1]) - 1
        return start <= end_node and st_node < end
    
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
            elif self.nodes_ref_list[i][2] >= end:
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
            start, end = map(int, (start, end))
            if len(seq) > 12 or is_telomere_motif(seq):
                continue
            if mem and start - mem[-1][2] <= 50 and contig == mem[-1][0]:
                p_contig, p_start, p_end, p_seq_l = mem.pop()
                mem.append((p_contig, p_start, end, p_seq_l + [seq]))
                continue
            mem.append((contig, start, end, [seq]))
        for contig, start, end, seq_l in mem:
            nodes_sequence, left_flanking, right_flanking = self.binary_search_nodes(start, end)

            #if (start, end) in [(10911966, 10911975), (11018316, 11018333)]:
                #print(start, end, '\n', nodes_sequence, '\n', '\n'.join(map(lambda x : self.nodes[x], nodes_sequence)))
            if len(nodes_sequence) > 1:
                self.repeat_nodes[contig, start, end] = nodes_sequence
                self.repeat_motif[contig, start, end] = seq_l
                left_split = start - self.nodes_ref_ind[left_flanking]
                right_split = end - self.nodes_ref_ind[right_flanking]
                self.flanking_nodes[contig, start, end] = (left_flanking, left_split, right_flanking, right_split)
            #else:
                #print(self.c)
        #print(self.repeat_nodes.keys(), self.database_ref.keys())
        mem.clear()
        file.close()

    def _entropy(self, counter: Counter) -> float:
        counts = np.array(list(counter.values()), dtype=float)
        probs = counts / counts.sum()
        entropy = -np.sum(probs * np.log2(probs))
        return float(entropy)

    def _move_left_flank(self, first_node_ind: int, walk: List[int], left_flank: str) -> Tuple[str, str]:
        moving_ind: int = 0
        i_beg: int = first_node_ind
        len_left: int = len(left_flank)
        parts_list: deque[str] = deque()
        while True:
            # print(moving_ind)
            while (len_left - moving_ind) < 50 and i_beg >= 0:
                len_left += len(self.nodes[abs(walk[i_beg])])
                #left_flank = self.rev_if_neg(walk[i_beg]) + left_flank
                parts_list.appendleft(self.rev_if_neg(walk[i_beg]))
                i_beg -= 1
            if i_beg < 0 and moving_ind >= len(left_flank):
                break
            left_flank = "".join(parts_list) + left_flank
            parts_list.clear()
            tuples_counts = Counter(
                [left_flank[i:i + 3] for i in range((len_left - moving_ind) - 10, (len_left - moving_ind) - 2)]
            )
            entropy = self._entropy(tuples_counts)
            if entropy < 2.77:
                moving_ind += 1
            else:
                break
        boundary = len_left - moving_ind
        left_flank_final = left_flank[boundary - 50:boundary]
        repeat_prefix = left_flank[boundary:]
        return left_flank_final, repeat_prefix

    def _move_right_flank(self, first_node_ind: int, walk: List[int], right_flank: str) -> Tuple[str, str]:
        moving_ind: int = 0
        i_end: int = first_node_ind
        len_right: int = len(right_flank)
        parts_list: deque[str] = deque()
        while True:
            # print(moving_ind)
            while (len_right - moving_ind) < 50 and i_end < len(walk):
                len_right += len(self.nodes[abs(walk[i_end])])
                parts_list.append(self.rev_if_neg(walk[i_end]))
                i_end += 1
            # print(right_flank)
            if i_end == len(walk) and moving_ind >= len(right_flank):
                break
            right_flank += "".join(parts_list)
            parts_list.clear()
            tuples_counts = Counter([right_flank[i:i + 3] for i in range(moving_ind, moving_ind + 8)])
            entropy = self._entropy(tuples_counts)
            if entropy < 2.77:
                moving_ind += 4
            else:
                break
        return right_flank[moving_ind:moving_ind + 50], right_flank[:moving_ind]

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
            genome: str,
            ) -> None:
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
        #if (start, end) in [(10911966, 10911975), (11018316, 11018333)]:
            #print(start, end, '\n', stack, '\n', '\n'.join(map(lambda x : self.nodes[x], stack)))
        node_spans: List[Tuple[int, int, int]] = []
        cursor: int = 0
        seq: str = ""
        for x in stack:
            s = self.rev_if_neg(x)
            seq += s
            node_spans.append((abs(x), cursor, cursor + len(s)))
            cursor += len(s)
        flanks: Tuple[int, int, int, int] = self.flanking_nodes[ref_key]
        left, repeat, right = self.split_by_flanks(seq, node_spans, flanks)
        #if (start, end) in [(10911966, 10911975), (11018316, 11018333)]:
            #print(start, end, '\n', flanks)
        #if current_start == 127453829:
            #print(f"start {current_start}")
        if i - 1 >= 0:
            new_left, rep_add = self._move_left_flank(i - 1, walk, left)
            left = new_left
            repeat = rep_add + repeat
            #print("left")
        if i + len(stack) < len(walk):
            new_right, rep_add = self._move_right_flank(i + len(stack), walk, right)
            repeat += rep_add
            right = new_right
            #print("right")
        current_start = current + len(left)
        current_end = current_start + len(repeat)
        if (ref_contig, start, end) not in paths:
            paths[ref_contig, start, end] = {genome: (current_start, current_end)}
            paths_seqs[ref_contig, start, end] = {genome: (left, repeat, right)}
        elif genome not in paths[ref_contig, start, end] or (paths[ref_contig, start, end][genome][1] -
                paths[ref_contig, start, end][genome][0]) < current_end - current_start:
            paths[ref_contig, start, end][genome] = (current_start, current_end)
            paths_seqs[ref_contig, start, end][genome] = (left, repeat, right)


    def analysis(
            self,
            mismatch_penalty: int,
            match_score: int,
            threshold: int
    ):
        self.rev_if_neg = lambda x: self.nodes[-x][::-1] if x < 0 else self.nodes[x]
        self.nodes_ref_ind.clear()
        self.nodes_ref_list.clear()
        node_occurrences = defaultdict(list)
        for (gen, contig), nodes_dict in self.walks.items():
            for (start_w, _), walk in nodes_dict.items():
                start_n = start_w
                for i, node in enumerate(walk):
                    node_occurrences[node].append((walk, i, start_n, gen))
                    start_n += len(self.nodes[abs(node)])
        paths: Dict[Tuple[str, int, int], Dict[str, Tuple[int, int]]] = {}
        paths_seqs: Dict[Tuple[str, int, int], Dict[str, Tuple[str, str, str]]] = {}
        for ref_key, nodes_rf in self.repeat_nodes.items():
            for node in nodes_rf:
                for walk, i, start_n, gen in node_occurrences.get(node, []):
                    self._evaluate(
                        i=i,
                        walk=walk,
                        threshold=threshold,
                        match_score=match_score,
                        mismatch_penalty=mismatch_penalty,
                        nodes_rf=nodes_rf,
                        current=start_n,
                        ref_key=ref_key,
                        paths=paths,
                        paths_seqs=paths_seqs,
                        genome=gen
                    )
        with open(self.l, 'w') as log:
            for key, info in paths_seqs.items():
                blocks: List[Tuple[str, str, str, str]] = []
                ref_ind: Optional[int] = None
                ind: int = 0
                for genome, (left, repeat, right) in info.items():
                    if genome == self.reference:
                        ref_ind = ind
                    blocks.append((left, repeat, right, genome))
                    ind += 1
                del ind
                if ref_ind is not None:
                    blocks[0], blocks[ref_ind] = blocks[ref_ind], blocks[0]
                print(
                    f"Tandem repeat motif 5'-{"-*-".join(self.repeat_motif[key])}"
                    f"-3'\nIn {key[0]} at [ {paths[key][self.reference][0]}, {paths[key][self.reference][1]} ]:",
                    file=log
                )
                self.pretty_print(blocks, key, paths, log=log)
            #print(len(paths))

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
        return left, repeat, right

    def pretty_print(self, blocks: List[Tuple[str, str, str, str]], ref_key: Tuple[str, int, int],
                     paths: Dict[Tuple[str, int, int], Dict[str, Tuple[int, int]]], log):
        max_rep: int = max(len(t) for _, t, _, _ in blocks)
        max_left: int = max(len(l) for l, _, _, _ in blocks)
        for left, rep, right, gen in blocks:
            approx_start: int = paths[ref_key][gen][0]
            approx_end: int = paths[ref_key][gen][1]
            if gen != self.reference:
                print(f"Found in {gen} at [ {approx_start}, {approx_end} ]:", file=log)
            print(f"{left:>{max_left}}  {rep:<{max_rep}}  {right}", file=log)
        print(file=log)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', "--repeats", required=True, help="tandem-repeats")
    parser.add_argument('-l', "--log", required=True, help="logfile")
    parser.add_argument('-c', "--chr", required=True, help="logfile")
    parser.add_argument("--ref", required=True, help="reference genome")
    args = parser.parse_args()
    r, l, c, ref = args.repeats, args.log, args.chr, args.ref
    processing = Processing(r, l, c, ref)
    #import time
    #start = time.time()
    processing.analysis(-1, 1, 2)
    #print("Computation ended in %f seconds" % (time.time() - start))
