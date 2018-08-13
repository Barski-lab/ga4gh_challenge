cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.sort_input,
                  "entryname": inputs.sort_input.basename,
                  "writable": true
                }
              ]
    }


- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() {
      if (inputs.csi && !inputs.bai){
        return '.csi';
      } else {
        return '.bai';
      }
    };
  - var default_bam = function() {
      if (inputs.trigger == true){
        return inputs.sort_input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".bam";
      } else {
        return inputs.sort_input.location.split('/').slice(-1)[0];
      }
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:

  bash_script_sort:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        echo "Run: samtools sort " ${@:1}
        samtools sort "${@:1}"
      else
        echo "Skip samtools sort " ${@:1}
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run samtools sort with all input parameters or skip it if trigger is false

  bash_script_index:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        echo "Run: samtools index " ${@:1}
        samtools index "${@:1}"
      else
        echo "Skip samtools index " ${@:1}
      fi
    inputBinding:
      position: 20
    doc: |
      Bash function to run samtools index with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    doc: |
      If true - run samtools, if false - return sort_input and optional index file in secondaryFiles, previously staged
      into output directory.

  sort_compression_level:
    type: int?
    inputBinding:
      position: 11
      prefix: -l
    doc: |
      SORT: desired compression level for the final output file, ranging from 0 (uncompressed)
      or 1 (fastest but minimal compression) to 9 (best compression but slowest to write),
      similarly to gzip(1)'s compression level setting.
      If -l is not used, the default compression level will apply.

  sort_output_filename:
    type: string?
    inputBinding:
      position: 12
      prefix: -o
      valueFrom: |
        ${
            if (self == "" || inputs.trigger == false){
              return default_bam();
            } else {
              return self;
            }
        }
    default: ""
    doc: |
      Write the final sorted output to FILE. Only out.bam|out.cram.
      If output file extension is set to SAM, tool will fail on the index step

  threads:
    type: int?
    doc: |
      Set number of sorting and compression threads [1] (Only for sorting)

  sort_input:
    type: File
    inputBinding:
      position: 16
    doc: |
      Input only in.sam|in.bam|in.cram. Optionally could be supplemented with index file in secondaryFiles

  csi_interval:
    type: int?
    inputBinding:
      position: 24
      prefix: -m
    doc: |
      Set minimum interval size for CSI indices to 2^INT [14]

  csi:
    type: boolean?
    doc: |
      Generate CSI-format index for BAM files. If input isn't cram.

  bai:
    type: boolean?
    doc: |
      Generate BAI-format index for BAM files [default]. If input isn't cram.

outputs:
  bam_bai_pair:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.sort_output_filename == "" || inputs.trigger == false){
              return default_bam();
            } else {
              return inputs.sort_output_filename;
            }
        }
    secondaryFiles:
      ${
          if (inputs.trigger == true){
            return self.basename + ext();
          } else {
            return inputs.sort_input.secondaryFiles?inputs.sort_input.secondaryFiles:"null";
          }
      }

baseCommand: [bash, '-c']

arguments:
#   run script sort position 5
  - valueFrom: |
      ${ return inputs.trigger ? "true" : "false" }
    position: 6
  # -l - position 11
  # -o sort_output_filename - position 12
  - valueFrom: bam
    position: 13
    prefix: -O
    # -n - position 14
  - valueFrom: $(inputs.threads?inputs.threads:1)
    position: 15
    prefix: -@
  # sort_input - position 16
  - valueFrom: ";"
    position: 17
    shellQuote: false

  - valueFrom: "bash"
    position: 18
  - valueFrom: "-c"
    position: 19
#   run script index position 20
  - valueFrom: |
      ${ return inputs.trigger ? "true" : "false" }
    position: 21
  - valueFrom: $(inputs.bai?'-b':inputs.csi?'-c':[])
    position: 23
    # -m - position 24
  - valueFrom: $(inputs.threads?inputs.threads:1)
    position: 25
    prefix: -@
  - valueFrom: |
      ${
          if (inputs.sort_output_filename == "" || inputs.trigger == false){
            return default_bam();
          } else {
            return inputs.sort_output_filename;
          }
      }
    position: 26
  - valueFrom: |
      ${
          if (inputs.sort_output_filename == "" || inputs.trigger == false){
            return default_bam() + ext();
          } else {
            return inputs.sort_output_filename + ext();
          }
      }
    position: 27
