#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'csv'
require 'bio'

$:.push File.join(File.dirname(__FILE__),'..','bioruby-iotu','lib')
require 'bio-iotu' #Has OtuTable class

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :img_basedir => '/srv/whitlam/bio/db/img/3.5/genomes/finished/',
    :out_fasta => 'grinder_reference.fna',
    :out_abundances => 'grinder_abundances.tsv',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <arguments>
      
      Convert the abundance and species information from mince.rb and output the grinder_reference.fna file and 
        \n"
      
    opts.on("-m", "--mince MINCE_OUTPUT", "Output from mince.rb [required]") do |arg|
      options[:mince_file] = arg
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0 or options[:mince_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  File.open(options[:out_fasta],'w') do |out_fasta|
    File.open(options[:out_abundances],'w') do |out_abundances|
      CSV.foreach(options[:mince_file], :col_sep => "\t") do |row|
        raise unless row.length > 1
        abundance = row[0].to_f
        img_identifier = row[1]
        raise unless img_identifier.match(/^\d+$/)
        
        # Concatenate the fna file from that OTU into the references.fna file, and output the abundances
        fasta = Bio::FlatFile.foreach(File.join(options[:img_basedir],img_identifier,"#{img_identifier}.fna")) do |seq|
          out_fasta.puts seq.entry
          out_abundances.puts [
            seq.definition.split(/\s/)[0],
            abundance,
          ].join(' ')
        end
      end
    end
  end
end #end if running as a script