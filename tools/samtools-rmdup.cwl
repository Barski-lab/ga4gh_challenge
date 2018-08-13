cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          return inputs.bam_file.location.split('/').slice(-1)[0];
        };
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.bam_file,
                  "entryname": "input_file_backup",
                  "writable": true
                }
              ]
    }

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4

inputs:

  bash_script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        samtools rmdup "${@:1}"
      else
        echo "Skip samtools rmdup " ${@:1}
        mv ${@: -2}
      fi
    inputBinding:
      position: 1
    doc: |
      Bash function to run samtools rmdup with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 2
      valueFrom: |
        ${ return self ? "true" : "false" }
    doc: |
      If true - run samtools rmdup, if false - return bam_file, previously staged into output directory and optional index file
      Use valueFrom to return string instead of boolean, because if return boolean False, argument is not printed

  single_end:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 3
      prefix: '-s'
    doc: |
      rmdup for SE reads

  force_single_end:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 4
      prefix: '-S'
    doc: |
      treat PE reads as SE in rmdup (force -s)

  bam_file:
    type: File
    inputBinding:
      position: 10
    doc: |
      Input sorted bam file (index file is optional)

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 11
      valueFrom: |
        ${
            if (self == "" || inputs.trigger == false){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: ""
    doc: |
      Writes the output bam file to output_filename if set,
      otherwise generates output_filename on the base of bam_file

outputs:
  rmdup_output:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == "" || inputs.trigger == false){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }
    secondaryFiles: |
      ${
          if (inputs.bam_file.secondaryFiles && inputs.trigger == false){
            return inputs.bam_file.secondaryFiles;
          } else {
            return "null";
          }
        }
    doc: File with removed duplicates or bam_file with optional secondaryFiles

  rmdup_log:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == "" || inputs.trigger == false){
            return default_output_filename().split('.').slice(0,-1).join('.') + '.rmdup';
          } else {
            return inputs.output_filename.split('.').slice(0,-1).join('.') + '.rmdup';
          }
        }

baseCommand: [bash, '-c']

arguments:
  - valueFrom: |
      ${
        if (inputs.output_filename == "" || inputs.trigger == false){
          return " > " + default_output_filename().split('.').slice(0,-1).join('.') + ".rmdup 2>&1";
        } else {
          return " > " + inputs.output_filename.split('.').slice(0,-1).join('.') + ".rmdup 2>&1";
        }
      }
    position: 100000
    shellQuote: false
