#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'reachable'
require 'pp'
 
$:.push File.join(File.dirname(__FILE__),'..','bioruby-iotu','lib')
require 'bio-iotu' #Has OtuTable class

$:.push File.join(File.dirname(__FILE__),'..','bioruby-taxonomy_definition_files','lib')
require 'bio-taxonomy_definition_files'

# Extend this class with some case-specific methods
class Bio::IMG::TaxonomyDefinitionFile
  def get_random_genus_with_at_least_two_species
    counts = {}
    each do |tax|
      counts[tax.genus] ||= 0
      counts[tax.genus] += 1
    end
    counts.select{|genus, count| count >= 2}.to_a.sample[0]
  end
  
  def get_two_random_orders
    counts = {}
    each do |tax|
      counts[tax.order] ||= 0
      counts[tax.order] += 1
    end
    counts.select{|genus, count| count >= 2}.keys.shuffle[0..1]
  end
end

class Array
  def low_medium_high_abundance
    medium = 0.05
    low = 0.01
    high = 0.25
    total_abundance = inject{|sum,i|sum+=i}
    return [
      total_abundance*low,
      total_abundance*medium,
      total_abundance*high,
    ].collect{|s| s.to_i}
  end
end

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :otu_table_columns => [0],
    :same_genus1 => false,
    :same_genus2 => false,
    :different_orders1 => false,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Take an OTU table and a list of sequenced genomes, output a list of abundances for each species.
        \n"
      
    opts.on("-o", "--otu-table OTU_TABLE_FILE", "OTU table to be simulated [required]") do |arg|
      options[:otu_table] = arg
    end
    opts.on("-t", "--possible-taxonomies POSSIBLE_TAXONOMIES_FILE", "Path to IMG taxonomy file i.e. IMG3.5_release_metadata.tsv [required]") do |arg|
      options[:sequenced_genomes_taxonomy] = arg
    end
    opts.on("-c", "--otu-table-columns COLUMNS", "Columns [default: #{options[:otu_table_columns].join(',')}]") do |arg|
      options[:otu_table_columns] = arg.split(",").collect{|s| s.strip.to_i}
    end
    opts.on("--same-genus-different-coverages1", "Take a random genus and include exactly 2 species, with coverage profiles up,middle,up vs. down,middle,down [default: #{options[:same_genus1]}]") do
      options[:same_genus1] = true
    end
    opts.on("--same-genus-different-coverages2", "Take a random genus and include exactly 2 species, with coverage profiles middle,middle,up vs. middle,middle,down [default: #{options[:same_genus2]}]") do
      options[:same_genus2] = true
    end
    opts.on("--different-order-same-coverages1", "Take a random pair of orders and include exactly 1 species from each order, and then give them each the same coverage profile [default: #{options[:different_orders1]}]") do
      options[:different_orders1] = true
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
  
  # Sorted in descending abundance in the OTU table
  table = Bio::OtuTable.read(options[:otu_table])
  log.debug "Read #{table.otu_identifiers.length} OTUs spread across #{table.samples.length} samples in the OTU table"
  
  # Read the list of sequenced genomes file
  log.info "Reading genomes to model from #{options[:sequenced_genomes_taxonomy]}"
  sequenced_genomes = Bio::IMG::TaxonomyDefinitionFile.read(options[:sequenced_genomes_taxonomy])
  
  # Remove from the table all except the samples to be processed
  log.info "Using abundances from samples #{options[:otu_table_columns].join(",")}"
  samples_to_keep_names = options[:otu_table_columns].collect{|i| table.samples.keys[i]}
  table.samples = table.samples.select do |name, abundances|
    samples_to_keep_names.include? name
  end
  table.remove_empty_rows!
  log.debug "After trimming the OTU table down to remove species not present, the OTU table has these samples: #{table.samples.keys.join(', ')}, and #{table.otu_identifiers.length} OTUs"
  
  # Add in the spike-ins. First need to know the median, Q1, Q3
  median_abundances = {}
  q1_abundances = {}
  q3_abundances = {} 
  table.samples.each do |sample, abundances|
    q1_abundances[sample], median_abundances[sample], q3_abundances[sample] = abundances.low_medium_high_abundance
    
    log.info "For sample '#{sample}', found low,medium,high abundances #{q1_abundances[sample]}, #{median_abundances[sample]}, #{q3_abundances[sample]}"
  end
  
  spiking_in = false
  if (options[:same_genus1] or options[:same_genus2] or options[:different_orders1])
    spiking_in = true
  end
  if spiking_in and table.samples.length != 3
    log.error "Need exactly 3 samples for the current spike-ins, I found #{table.samples.length}"
  end
  
  spike_ins = {}# IMG-definition => {sample => abundances}
  sample_names = table.sample_names
  sample_names.each do |name|
    spike_ins[name] = {}
  end
  
  # Take a random genus and include exactly 2 species, with coverage profiles up,middle,up vs. down,middle,down
  if options[:same_genus1]
    # Pick a random genus
    random_genus = sequenced_genomes.get_random_genus_with_at_least_two_species
    possible_species_to_spike_in = sequenced_genomes.select{|tax| tax.genus == random_genus}
    possible_species_to_spike_in.shuffle!
    random_species1 = possible_species_to_spike_in[0]
    random_species2 = possible_species_to_spike_in[1]
    
    log.info "Spiking in #{random_species1.genus_species} and #{random_species2.genus_species} for spike in #1"
    
    # Assign the spiked-in abundances
    spike_ins[sample_names[0]][random_species1] = q3_abundances[sample_names[0]]
    spike_ins[sample_names[1]][random_species1] = median_abundances[sample_names[1]]
    spike_ins[sample_names[2]][random_species1] = q3_abundances[sample_names[2]]
    
    spike_ins[sample_names[0]][random_species2] = q1_abundances[sample_names[0]]
    spike_ins[sample_names[1]][random_species2] = median_abundances[sample_names[1]]
    spike_ins[sample_names[2]][random_species2] = q1_abundances[sample_names[2]]
    
    # Remove anything from that genus to the abundances file
    sequenced_genomes.reject! do |tax|
      tax.genus == random_genus
    end
  end
  log.info "After spike-in #1, there are #{sequenced_genomes.length} different OTUs available for random modelling"
  
  # Take a random genus and include exactly 2 species, with coverage profiles middle,middle,up vs. middle,middle,down
  if options[:same_genus2]
    # Pick a random genus
    random_genus = sequenced_genomes.get_random_genus_with_at_least_two_species
    possible_species_to_spike_in = sequenced_genomes.select{|tax| tax.genus == random_genus}
    possible_species_to_spike_in.shuffle!
    random_species1 = possible_species_to_spike_in[0]
    random_species2 = possible_species_to_spike_in[1]
    
    log.info "Spiking in #{random_species1.genus_species} and #{random_species2.genus_species} for spike in #2"
    
    # Assign the spiked-in abundances
    spike_ins[sample_names[0]][random_species1] = median_abundances[sample_names[0]]
    spike_ins[sample_names[1]][random_species1] = median_abundances[sample_names[1]]
    spike_ins[sample_names[2]][random_species1] = q3_abundances[table.samples[2]]
    
    spike_ins[sample_names[0]][random_species2] = median_abundances[sample_names[0]]
    spike_ins[sample_names[1]][random_species2] = median_abundances[sample_names[1]]
    spike_ins[sample_names[2]][random_species2] = q1_abundances[sample_names[2]]
    
    # Remove anything from that genus to the abundances file
    sequenced_genomes.reject! do |tax|
      tax.genus == random_genus
    end
  end
  log.info "After spike-in #2, there are #{sequenced_genomes.length} different OTUs available for random modelling"
  
  if options[:different_orders1]
    random_order1, random_order2 = sequenced_genomes.get_two_random_orders
    log.info "Chose random orders #{random_order1} and #{random_order2} for options[:different_orders]"
    
    log.info "Spiking in #{random_species1.genus_species} and #{random_species2.genus_species} for spike in #3"
    
    # Assign the spiked-in abundances
    spike_ins[sample_names[0]][random_species1] = median_abundances[sample_names[0]]
    spike_ins[sample_names[1]][random_species1] = median_abundances[sample_names[1]]
    spike_ins[sample_names[2]][random_species1] = median_abundances[sample_names[2]]
    
    spike_ins[sample_names[0]][random_species2] = median_abundances[sample_names[0]]
    spike_ins[sample_names[1]][random_species2] = median_abundances[sample_names[1]]
    spike_ins[sample_names[2]][random_species2] = median_abundances[sample_names[2]]
    
    # Remove anything from that genus to the abundances file
    sequenced_genomes.reject! do |tax|
      tax.order == random_order1 or tax.order == random_order2
    end
  end
  log.info "After the spike-ins, there are #{sequenced_genomes.length} different OTUs available for random modelling"
  
  pp spike_ins
  
  # Assign the rest of the abundances to each of the OTU ids available in the table
  sequenced_genomes.shuffle!
  if table.otu_identifiers.length > sequenced_genomes.length
    raise "There are more OTUs in the OTU table (#{table.otu_identifiers.length}) than genomes available for modelling (#{sequenced_genomes.length}), so quitting. Get more!"
  end
  assigned_genomes = sequenced_genomes[0...table.otu_identifiers.length]
  
  table.samples.each do |sample, sample_abundances|
    abundances_filename = "abundances.#{sample}.csv"
    log.info "Now writing abundances from sample #{sample} to #{abundances_filename}"
    
    File.open(abundances_filename,'w') do |f|
      # Write in the spike-in data
      spike_ins[sample].each do |taxon, abundance|
        f.puts [
          abundance,
          taxon.definition_line
        ].join("\t")
      end
      
      # Write the regular abundances
      sample_abundances.each_with_index do |abundance, i|
        f.puts [
          abundance,
          assigned_genomes[i].definition_line,
        ].join("\t")
      end
    end
  end
  
end #end if running as a script