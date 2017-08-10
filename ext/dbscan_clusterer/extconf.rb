require 'mkmf'

$LDFLAGS += " -lm"

extension_name = 'dbscan_clusterer'
dir_config(extension_name)
create_makefile(extension_name)
