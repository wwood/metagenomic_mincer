#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'

$:.push File.join(File.dirname(__FILE__),'..','bioruby-taxonomy_definition_files','lib')
require 'bio-taxonomy_definition_files' #has IMG taxonomy parser file

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Takes in an IMG taxonomy file, and removes from it all but one strain from each species\n\n"
      
    opts.on("-t", "--img-taxonomies-file PATH", "Path to IMG taxonomy file i.e. IMG3.5_release_metadata.tsv [required].") do |arg|
      options[:full_img_taxonomies_file] = arg
    end
    
    # logger options
    opts.separator "\nVerbosity:\n\n"
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end; o.parse!
  if ARGV.length != 0
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  # Read in the IMG taxonomy
  
  taxonomies = Bio::IMG::TaxonomyDefinitionFile.read(options[:full_img_taxonomies_file])
  raise unless taxonomies.length > 1000
  
  # Randomise so that there is no biases in the order of this file. That's unlikely, but still, good to do.
  taxonomies.shuffle!
  
  log.info "Read in #{taxonomies.length} species to work with"
  
  accepted_genus_species = []
  taxonomies.reject! do |lineage|
    index = "#{lineage.genus}_#{lineage.species}"
    if accepted_genus_species.include?(index)
      true
    else
      accepted_genus_species.push index
      false
    end
  end
  
  log.info "After removing different strains of the same species, #{taxonomies.length} different species were left"
  
  # Headers
  puts %w(taxon_oid
Domain
Status
Genome Name
Phylum
Class
Order
Family
Genus
Species
Strain
Release Date
IMG Release
).join("\t")
  # Data
  taxonomies.each do |lineage|
    puts lineage.definition_line
  end
end #end if running as a script