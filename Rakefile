require 'rake/extensiontask'


# somewhere in your Rakefile, define your gem spec
spec = Gem::Specification.new do |s|
  s.name    = "dbscan_clusterer"
  s.version = "0.1.0"
  s.summary = "dbscan for Ruby"
  s.author  = "Julien Boulnois"

  s.platform = Gem::Platform::RUBY
  s.extensions = FileList["ext/**/extconf.rb"]

  s.files = Dir.glob("ext/**/*.{c,rb}") +
            Dir.glob("lib/*.{rb,so}")
end

Gem::PackageTask.new(spec) do |pkg|
end

# feed the ExtensionTask with your spec
Rake::ExtensionTask.new('dbscan_clusterer', spec)