require 'rake/extensiontask'

spec = Gem::Specification.load("dbscan_clusterer.gemspec")

task :default => [:compile]

Gem::PackageTask.new(spec) do |pkg|
end

# feed the ExtensionTask with your spec
Rake::ExtensionTask.new('dbscan_clusterer', spec)