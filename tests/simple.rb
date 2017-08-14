require 'dbscan_clusterer'
require 'pp'

data = [[ 0, 10 ],
        [ 0, 11 ],
        [ 0, 12 ],
        [ 20, 33 ],
        [ 21, 32 ],
        [ 59, 77 ],
        [ 58, 79 ],
        [ 58, 76 ],
        [ 300, 70 ],
        [ 500, 300 ],
        [ 500, 302 ]]

pp DbscanClusterer.dbscan(data,4,1)
pp DbscanClusterer.dbscan(data,4,1,:euclidean2d)
pp DbscanClusterer.dbscan(data,4,1,:ruby,
        Proc.new{|a,b| Math.sqrt(((a[0] - b[0]) ** 2) + ((a[1] - b[1]) ** 2))})