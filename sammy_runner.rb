#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'bio'
require 'progressbar'
require 'parallel'
require 'reachable'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :img_basedir => '/srv/whitlam/bio/db/img/3.5/genomes/finished/',
    :out_fasta => 'grinder_reference.fna',
    :out_abundances => 'grinder_abundances.tsv',
    :threads => 1,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Convert the abundance and species information from mince.rb and output the grinder_reference.fna file and 
        \n"
      
    opts.on("-m", "--mince MINCE_OUTPUT", "Output from mince.rb [required]") do |arg|
      options[:mince_file] = arg
    end
    opts.on("--img-basedir DIRECTORY", "Base directory to IMG fasta files [default: #{options[:img_basedir]}]") do |arg|
      options[:img_basedir] = arg
    end
    opts.on("-n","--number READS", "Number of reads to simulate in total [required]") do |arg|
      options[:total_number_of_reads] = arg
    end
    opts.on("-t","--num-threads NUM_THREADS", "Number of reads to simulate in total [required]") do |arg|
      options[:threads] = arg.to_i
      raise unless options[:threads] > 0
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0 or options[:mince_file].nil? or options[:img_basedir].nil? or options[:total_number_of_reads].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  # First work out the total number of OTUs attributed to each
  total_otu_count = 0
  otus = []
  class Otu
    attr_accessor :abundance, :img_identifier
  end
  CSV.foreach(options[:mince_file], :col_sep => "\t") do |row|
    total_otu_count += row[0].to_i

    otu = Otu.new
    otu.abundance = row[0].to_f
    otu.img_identifier = row[1]
    otus.push otu
  end
  non_zero_otus = otus.select{|a| a.abundance>0}
  log.info "Total OTU count was #{otus.length}, total with non-zero abundance count #{non_zero_otus.length}"
  log.info "Total abundance was #{otus.reach.abundance.reduce(:+)}"
  
  progress = ProgressBar.new('collect', `wc -l #{options[:mince_file]}`.to_i)
  Parallel.each_with_index(non_zero_otus, :in_threads => options[:threads]) do |otu, index|
    fasta_file = File.join(options[:img_basedir],otu.img_identifier,"#{otu.img_identifier}.fna")
    number_of_reads = (options[:total_number_of_reads].to_f*otu.abundance/total_otu_count).round
    log.debug "From abundance #{otu.abundance}, the total number of reads is #{number_of_reads}"
    cmd = [
      'sammy.pl',
      '-r',
      fasta_file,
      '-n',
      number_of_reads,
      ">otu#{index}.fa"
    ].join(' ')
    `#{cmd}`
    progress.inc
  end
  progress.finish
end #end if running as a script