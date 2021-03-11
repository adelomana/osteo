 '''
    # plot zscores as a matrix
    matplotlib.pyplot.imshow(zscores_regulon, aspect='auto')
    matplotlib.pyplot.colorbar()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.show()

    # run some clustering
    #clustering = OPTICS(min_samples=5).fit_predict(zscores_regulon)
    #print(clustering, len(clustering))

    clustering = AgglomerativeClustering().fit_predict(zscores_regulon)
    print(clustering)

    for i in range(len(zscores_regulon)): 
        if clustering[i] == 0:
            matplotlib.pyplot.plot(time_trajectory, zscores_regulon[i], 'o-', alpha=0.5)
        matplotlib.pyplot.xlabel('Time (day)')
        matplotlib.pyplot.ylabel('zscore')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.show()

    '''







pca = sklearn.decomposition.PCA(n_components=5)
    projection = pca.fit_transform(numpy.transpose(zscores_regulon))

    PC1_loadings = pca.components_[0]

    PC1_loadings.sort()
    print('loadings', PC1_loadings, len(PC1_loadings))

    for i in range(len(PC1_loadings)):
        if PC1_loadings[i] > 0:
            the_color = 'red'
        else:
            the_color = 'blue'
        matplotlib.pyplot.plot(i, numpy.abs(PC1_loadings[i]), 'o-')
    matplotlib.pyplot.show()

    ### plot PCA
    for i in range(len(projection)):
        matplotlib.pyplot.plot(projection[i][0], projection[i][1], 'o', color=theColors[i], ms=20, alpha=0.5)
        matplotlib.pyplot.text(projection[i][0], projection[i][1], str(time_trajectory[i]), horizontalalignment='center', verticalalignment='center')
    matplotlib.pyplot.xlabel('PC1 ({}%)'.format(int(pca.explained_variance_ratio_[0]*100)))
    matplotlib.pyplot.ylabel('PC2 ({}%)'.format(int(pca.explained_variance_ratio_[1]*100)))
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.show()
