# somewhere in your Rakefile, define your gem spec
spec = Gem::Specification.new do |s|
  s.name    = "dbscan_clusterer"
  s.version = "0.2.0"
  s.summary = "dbscan for Ruby"
  s.author  = "Julien Boulnois"

  s.extensions = "ext/dbscan_clusterer/extconf.rb"

  s.files = Dir.glob("ext/dbscan_clusterer/*.{c,rb}") +
            Dir.glob("lib/*.{rb,so}")
end
