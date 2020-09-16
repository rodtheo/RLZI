import pandas as pd
import seaborn as sns
import sys
from pathlib import Path
import numpy as np

def max_subarray(numbers):
    """Find a contiguous subarray with the largest sum."""
    best_sum = 0  # or: float('-inf')
    best_start = best_end = 0  # or: None
    current_sum = 0
    for current_end, x in enumerate(numbers):
        if current_sum <= 0:
            # Start a new sequence at the current element
            current_start = current_end
            current_sum = x
        else:
            # Extend the existing sequence with the current element
            current_sum += x

        if current_sum > best_sum:
            best_sum = current_sum
            best_start = current_start
            best_end = current_end + 1  # the +1 is to make 'best_end' exclusive

    return best_sum, best_start, best_end

def q1(x):
    return x.quantile(0.25)

def q2(x):
    return x.quantile(0.75)

if __name__ == '__main__':

    dfs = []

    # print(sys.argv[1:])
    for f in sys.argv[1:]:
        fp = Path(f)
    # for i in range(1,24):

        # df = pd.read_table('/Users/rodtheo/Bioinfo/PROJETOS/STUDIES/LZ/Analises_paper_brazilian_bioinfo/asp/Ref_IFO6365chr{}/RLZ_k10_m10/synteny_to_draw.txt'.format(i), sep="\t")
        df = pd.read_table(f, sep="\t")

        df['strain'] = df.shape[0]*[fp.parent]

        dfs.append(df)

    df_all = pd.concat(dfs)

    df_all['len'] = df_all['End_2'] - df_all['Start_2']

    # ff = {'len': ['count']}
    # slstat = df_all.groupby('len').agg(ff)
    # print(slstat)
    # ##### uncomment to plot
    # # sns_plot = sns.distplot(df_all.loc[df_all['len'] < 16, 'len'], hist = True, kde = True, kde_kws = {'linewidth': 3})
    # sns_plot = sns.distplot(df_all.loc[df_all['len'] > 1, 'len'], hist = False, kde = True, kde_kws = {'linewidth': 3})
    # #
    #
    # # sns_plot.set_xlim(0, 30)
    # sns_plot.figure.savefig('density_synteny.png')
    #
    # sl = df_all.loc[df_all['len'] < 16, ]
    # ff = {'len': ['count']}
    # slstat = sl.groupby('len').agg(ff)
    # slstat['p'] = slstat/sl.shape[0]
    # print(slstat)
    # print(np.sum(slstat.iloc[:, 0]))
    # print(df_all.shape[0])
    # print(df_all.shape[0]-np.sum(slstat.iloc[:, 0]))
    #####

    f = {'len': ['median', 'std', q1, q2, 'max', 'count']}

    dfstat = df_all.groupby('strain').agg(f)
    df_all.sort_values(by=['len'], inplace=True, ascending=False)

    print(dfstat)

    subdfa = df_all.loc[df_all['fill'] == 'FF5733', :]
    subdf = subdfa.copy()
    subdf.reset_index(inplace=True)

    subdf.sort_values(by=['len'], inplace=True, ascending=False)

    print("==================================================")
    print("Top-10 Factors ordered by length")
    print(df_all.iloc[:10,:])
    print("==================================================")
    print("Top-10 Inverted Factors ordered by length")
    print(subdf.iloc[:10,:])
    print("==================================================")
    subdf.loc[subdf['len'] < 100, 'len'] = subdf.loc[subdf['len'] < 100, 'len'] * -1
    # # print(subdf)
    for f in sys.argv[1:]:
        fp = Path(f)
        subdfchr1 = subdf.loc[subdf['strain'] == fp.parent,:]
        # print(subdfchr1)
        subdfchr1ok = subdfchr1.sort_values(by=['Start_2'])
        # print(subdfchr1ok['len'].values)
        best_score, start, end = max_subarray(subdfchr1['len'].values)
        print("STRAIN=", str(fp.parent))
        print(best_score, start, end)
        print("==================================================")
        # print(subdfchr1.iloc[(start):(end), :])
    #.sort_values(by=['Start_2'])
