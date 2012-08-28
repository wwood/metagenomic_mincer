#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'reachable'
require 'csv'

$:.push File.join(File.dirname(__FILE__),'..','bioruby-iotu','lib')
require 'bio-iotu' #Has OtuTable class

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :otu_table_columns => [0],
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Take an OTU table and a list of sequenced genomes, output a list of abundances for each species.
        \n"
      
    opts.on("-o", "--otu-table OTU_TABLE_FILE", "OTU table to be simulated [required]") do |arg|
      options[:otu_table] = arg
    end
    opts.on("-t", "--possible-taxonomies POSSIBLE_TAXONOMIES_FILE", "list of species that are available to be picked [required]") do |arg|
      options[:sequenced_genomes_taxonomy] = arg
    end
    opts.on("-c", "--otu-table-columns COLUMNS", "Columns [default: #{options[:otu_table_columns].join(',')}]") do |arg|
      options[:otu_table_columns] = arg.split(",").collect{|s| s.strip.to_i}
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0 or options[:otu_table].nil? or options[:sequenced_genomes_taxonomy].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  # Read the list of sequenced genomes file
  sequenced_genomes = File.open(options[:sequenced_genomes_taxonomy]).readlines
  sequenced_genomes = sequenced_genomes.reject{|s| s.nil?}.reach.strip.to_a
  
  # Sorted in descending abundance in the OTU table
  table = Bio::OtuTable.read(options[:otu_table])
  log.debug "Read #{table.otu_identifiers.length} OTUs spread across #{table.samples.length} samples in the OTU table"
  
  # Read the list of sequenced genomes file
  sequenced_genomes = File.open(options[:sequenced_genomes_taxonomy]).readlines
  sequenced_genomes = sequenced_genomes.compact.reach.strip.to_a
  
  # Remove from the table all except the samples to be processed
  log.info "Using abundances from samples #{options[:otu_table_columns].join(",")}"
  samples_to_keep_names = options[:otu_table_columns].collect{|i| table.samples.keys[i]}
  table.samples = table.samples.select do |name, abundances|
    samples_to_keep_names.include? name
  end
  table.remove_empty_rows!
  log.debug "After trimming the OTU table down, the OTU table has these samples: #{table.samples.keys.join(', ')}, and #{table.otu_identifiers.length} OTUs"
  
  # Assign abundances to each of the OTU ids available in the table
  sequenced_genomes.shuffle!
  if table.otu_identifiers.length > sequenced_genomes.length
    raise "There are more OTUs in the OTU table (#{table.otu_identifiers.length}) than genomes available for modelling (#{sequenced_genomes.length}), so quitting. Get more!"
  end
  assigned_genomes = sequenced_genomes[0...table.otu_identifiers.length]
  
  table.samples.each do |sample, sample_abundances|
    log.info "Now writing abundances from sample #{sample}"
    
    File.open("abundances.#{sample}.csv",'w') do |f|
      sample_abundances.each_with_index do |abundance, i|
        f.puts [
          abundance,
          assigned_genomes[i],
        ].join("\t")
      end
    end
  end
  
end #end if running as a script