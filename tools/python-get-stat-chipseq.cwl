cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() {
        if (inputs.output_filename == null){
          let root = inputs.bowtie_log.basename.split('.').slice(0,-1).join('.');
          let ext = ".stat";
          return (root == "")?inputs.bowtie_log.basename+ext:root+ext;
        } else {
          return inputs.output_filename;
        }
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2


inputs:

  script:
    type: string?
    default: |
      #!/usr/bin/env python
      import sys, re
      TOTAL, ALIGNED, SUPRESSED, USED = 100, 80, 0, 0
      with open(sys.argv[1], 'r') as bowtie_log:
        for line in bowtie_log:
          if 'processed:' in line:
            TOTAL = int(line.split('processed:')[1])
          if 'alignment:' in line:
            ALIGNED = int(line.split('alignment:')[1].split()[0])
          if 'due to -m:' in line:
            SUPRESSED = int(line.split('due to -m:')[1].split()[0])
      USED = ALIGNED
      with open(sys.argv[2], 'r') as rmdup_log:
        for line in rmdup_log:
          if '/' in line and 'Skip' not in line:
            splt = line.split('/')
            USED = int((splt[1].split('='))[0].strip()) - int((splt[0].split(']'))[1].strip())
      print TOTAL, ALIGNED, SUPRESSED, USED
    inputBinding:
      position: 5
    doc: |
      Python script to get TOTAL, ALIGNED, SUPRESSED, USED values from log files

  bowtie_log:
    type: File
    inputBinding:
      position: 6
    doc: |
      Log file from Bowtie

  rmdup_log:
    type: File
    inputBinding:
      position: 7
    doc: |
      Log file from samtools rmdup

  output_filename:
    type:
    - "null"
    - string
    doc: |
      Name for generated output file


outputs:

  output_file:
    type: File
    outputBinding:
      glob: $(get_output_filename())

  total_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[0]))

  mapped_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[1]))

  supressed_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[2]))

  used_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_filename())
      outputEval: $(parseInt(self[0].contents.split(' ')[3]))

baseCommand: [python, '-c']
arguments:
  - valueFrom: $(" > " + get_output_filename())
    position: 100000
    shellQuote: false
