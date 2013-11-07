#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME}
      
      Concatenate and rename (by ordering) all the fasta files in the current directory, given that the fasta files are paired-end reads from sammy.\n\n"

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
  
  
  fastas_to_concat = Dir.glob('*.fa').sort
  log.info "Found #{fastas_to_concat.length} fasta files to concatenate"
  count = 1
  fastas_to_concat.each do |otu_file|
    Bio::FlatFile.foreach(otu_file) do |s|
      if count % 2 == 1
        puts ">#{(count+1)/2}_1"
        puts s.seq
      else
        puts ">#{count/2}_2"
        puts s.seq
      end
      count += 1
    end
  end
end #end if running as a script