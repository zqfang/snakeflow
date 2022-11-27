import seaborn as sns
aaas = ['red','green',...]
boxs = sns.boxplot(x='group',y='cellpct', data=cellpct2, palette=aaas, order=order,
                    showfliers=False, ax=ax[i])
# Add in points to show each observation    
# iterate over boxes
for m, box in enumerate(boxs.artists): // boxs.patches if matplotlib >=3.5
    box.set_edgecolor('black')
    box.set_facecolor('none')
    # iterate over whiskers and median lines
    for j in range(5*m,5*(m+1)):
        boxs.lines[j].set_color(aaas[m])