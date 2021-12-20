import argparse
import sys
import gc

import RNA

def hairpin_energy(fc, i, j):
    energy = fc.eval_hp_loop(i, j)
    
    return energy

def interior_energy(fc, i1, j1, i2, j2):
    energy = fc.eval_int_loop(i1, j1, i2, j2)

    return energy

def backtrack(Vbt, Wbt, dbn, i, j, arr):
    if arr == 0:    # array W
        if Wbt[i][j] is None:
            return

        if Wbt[i][j] == 'START':
            backtrack(Vbt, Wbt, dbn, i+1, j, 0)
        elif Wbt[i][j] == 'END':
            backtrack(Vbt, Wbt, dbn, i, j-1, 0)
        elif Wbt[i][j] == 'V':
            backtrack(Vbt, Wbt, dbn, i, j, 1)
        else:
            point = Wbt[i][j]
            backtrack(Vbt, Wbt, dbn, i, point, 0)
            backtrack(Vbt, Wbt, dbn, point+1, j, 0)
    else:           # array V
        if Vbt[i][j] is None:
            return

        dbn[i] = '('
        dbn[j] = ')'
        if Vbt[i][j] == 'HAIRPIN':
            return
        elif isinstance(Vbt[i][j], int):
            point = Vbt[i][j]
            backtrack(Vbt, Wbt, dbn, i+1, point, 0)
            backtrack(Vbt, Wbt, dbn, point+1, j-1, 0)
        else:
            i2, j2 = Vbt[i][j]
            backtrack(Vbt, Wbt, dbn, i2, j2, 1)

def zuker(seq):
    possible_pairs = {
        'G': 'C',
        'A': 'U',
        'U': 'A',
        'C': 'G'
    }

    n = len(seq)

    # dp array initialization
    V = [[0 for _ in range(n)] for _ in range(n)]
    W = [[0 for _ in range(n)] for _ in range(n)]
    Vbt = [[None for _ in range(n)] for _ in range(n)]
    Wbt = [[None for _ in range(n)] for _ in range(n)]
    dbn = ['.' for _ in range(n)]

    fc = RNA.fold_compound(seq)

    # base cases (pentanucleotide sequences)
    for i in range(n-4):
        j = i+4
        if possible_pairs[seq[i]] == seq[j]:
            V[i][j] = hairpin_energy(fc, i, j)
            Vbt[i][j] = 'HAIRPIN'
        else:
            V[i][j] = sys.maxsize

    for seqlen in range(5, n):
        for i in range(n-seqlen):
            gc.collect()
            j = i+seqlen
            
            # adding values of V
            if possible_pairs[seq[i]] == seq[j]:
                V_hairpin = hairpin_energy(fc, i, j)
                
                V_interior = None
                interior_point = None
                for i2 in range(i+1, j-1):
                    for j2 in range(i2+1, j):
                        if possible_pairs[seq[i2]] == seq[j2]:
                            val = interior_energy(fc, i, j, i2, j2) + V[i2][j2]
                            if V_interior is None or val < V_interior:
                                V_interior = val
                                interior_point = (i2, j2)
                
                if V_interior is None:
                    V_interior = sys.maxsize
                
                V_bifurcation = None
                bif_point = None
                for i2 in range(i+2, j-1):
                    val = W[i+1][i2] + W[i2+1][j-1]
                    if V_bifurcation is None or val < V_bifurcation:
                        V_bifurcation = val
                        bif_point = i2

                min_V = V[i][j] = min(V_hairpin, V_interior, V_bifurcation)
                if min_V == V_hairpin:
                    Vbt[i][j] = 'HAIRPIN'
                elif min_V == V_bifurcation:
                    Vbt[i][j] = bif_point
                else:
                    Vbt[i][j] = interior_point
            else:
                V[i][j] = sys.maxsize

            W_separate = None
            diff_point = None
            for i2 in range(i+1, j-1):
                val = W[i][i2] + W[i2+1][j]
                if W_separate is None or val < W_separate:
                    W_separate = val
                    diff_point = i2

            start_dangling = W[i+1][j]
            end_dangling = W[i][j-1]
            paired = V[i][j]
            min_W = W[i][j] = min(start_dangling, end_dangling, paired, W_separate)

            if min_W == start_dangling:
                Wbt[i][j] = 'START'
            elif min_W == end_dangling:
                Wbt[i][j] = 'END'
            elif min_W == paired:
                Wbt[i][j] = 'V'
            else:
                Wbt[i][j] = diff_point

    backtrack(Vbt, Wbt, dbn, 0, n-1, 0)
    dbnStr = ''.join(dbn)
    mfe = W[0][n-1] / 100.0

    return dbnStr, mfe

if __name__ == "__main__":
    RNA.params_load_RNA_Turner2004()
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq', dest='seq', type=str, required=True)

    args = parser.parse_args()

    db, mfe = zuker(args.seq)
    vDbn, vMfe = RNA.fold(args.seq)
    
    print('Format:\nOriginal Sequence\nZuker Structure\nViennaRNA structure\nMFE of both structures\n\n\n')
    print(f'{args.seq}\n{db}\n{vDbn}\nZuker MFE: {mfe} kcal/mol\t ViennaRNA MFE: {vMfe} kcal/mol')
