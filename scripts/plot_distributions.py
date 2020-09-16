import pandas as pd
import seaborn as sns
import sys
from pathlib import Path
import re
import matplotlib.pyplot as plt
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

if __name__=='__main__':

    plt.style.use('/Users/rodtheo/Bioinfo/PROJETOS/STUDIES/LZ/Analises_paper_brazilian_bioinfo/scripts/science.mplstyle')

    dfs = []

    # print(sys.argv[2:])
    ms = []
    for f in sys.argv[2:]:

        print(f)
        fpath = Path(f)

        mat = re.search(r"RLZ_k(\d+)_m(\d+)", f)
        print(mat.group(2))
        df = pd.read_table(f, sep="\t")
        # df['chr'] = df.shape[0]*['chr{}'.format(i)]
        ms.append(mat.group(2))
        df['m'] = df.shape[0]*[mat.group(2)]

        dfs.append(df)

    df_all = pd.concat(dfs)

    # print(df_all)

    df_all['len'] = df_all['End_2'] - df_all['Start_2']
    fig, ax = plt.subplots()
    # Turns off grid on the left Axis.
    ax.grid(False)
    # Turns off grid on the secondary (right) Axis.
    # ax.right_ax(False)
    df_all['loglen'] = np.log10(df_all['len'])
    for m in ms:

        subset = df_all[df_all['m'] == m]
        sns.distplot(np.log10(subset['len']), hist = False, kde = True,
                 kde_kws = {'linewidth': 3}, label=m, ax=ax)

    fig.canvas.draw()
    locs, labels = plt.xticks()
    # ax.set(xticklabels=[10 ** int(i.get_text().replace(u'\u2212', '-'))
                    # for i in labels])

    # sns_plot.figure.savefig(sys.argv[1])
    # Or for scientific notation:
    ax.set(xticklabels=["$10^" + i.get_text() + "$" for i in labels])
    ax.set_ylabel('Density')
    ax.set_xlabel('Phrase length')
    # with plt.style.context('science'):
    plt.savefig(sys.argv[1])

    fig, ax = plt.subplots()

    ax = sns.violinplot(x="m", y="loglen", data=df_all, order=[str(x) for x in sorted([int(x) for x in df_all['m'].unique()])])
    plt.savefig(sys.argv[1]+"violing.png")

    # re.match(r"RLZ_k(\d+)_m(\d+)", )

    # subdfa = df_all.loc[df_all['fill'] == 'FF5733', :]
    # subdf = subdfa.copy()
    # subdf.reset_index(inplace=True)
    #
    # # subdf.sort_values(by=['len'], inplace=True, ascending=False)
    #
    # subdf.loc[subdf['len'] < 70, 'len'] = subdf.loc[subdf['len'] < 70, 'len'] * -1
    # # print(subdf)
    # for i in range(1,24):
    #     subdfchr1 = subdf.loc[subdf['chr'] == 'chr{}'.format(str(i)),:]
    #     # print(subdfchr1)
    #     subdfchr1ok = subdfchr1.sort_values(by=['Start_2'])
    #     # print(subdfchr1ok['len'].values)
    #     best_score, start, end = max_subarray(subdfchr1['len'].values)
    #     print("chr", i)
    #     print(best_score, start, end)
    #     print("==================================================")
        # print(subdfchr1.iloc[(start):(end), :])
    #.sort_values(by=['Start_2'])
