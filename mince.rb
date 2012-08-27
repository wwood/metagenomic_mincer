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
  
  log.info "Using abundances from samples #{options[:otu_table_columns].join(",")}"
  options[:otu_table_columns].each do |sample_column_index|
    sample = table.samples.keys[sample_column_index] #just take the first one
    log.info "Using abundances from column #{sample_column_index} (#{sample})"
    
    all_abundances = table.samples[sample]
    sample_abundances = all_abundances.collect{|s| s.to_f}.reject{|s| s==0}.sort{|a,b| b<=>a}
    
    if sample_abundances.length > sequenced_genomes.length
      log.warn "There are more samples than possibilties for sequencing, so truncating the list of output genomes"
      sample_abundances = sample_abundances[0...sequenced_genomes.length]
    end
    log.info "Outputing abundances for #{sample_abundances.length} different OTUs from this sample"
    
    # Print the ID, and abundance
    sequenced_genomes.shuffle! #randomising this array means
    File.open("abundances.#{sample}.csv",'w') do |f|
      sample_abundances.each_with_index do |abundance, i|
        f.puts [
          abundance,
          sequenced_genomes[i],
        ].join("\t")
      end
    end
  end
  
end #end if running as a script