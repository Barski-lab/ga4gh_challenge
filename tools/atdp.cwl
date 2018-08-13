cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function(ext) {
        let root = inputs.input_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.input_file.basename+ext:root+ext;
    };
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.annotation_filename,
                  "entryname": inputs.annotation_filename.basename,
                  "writable": true
                }
              ]
    }


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/atdp:v0.0.1

inputs:

  script:
    type: string?
    default: |
        #!/bin/bash
        set -- "$0" "$@"
        echo "Original arguments:"
        for i in "$@";
            do echo $i;
        done;
        set -- "$1" --a=$(basename "${2:4}") "${@:3}"
        echo Updated arguments:
        for i in "$@";
            do echo $i;
        done;
        refgene-sort -i "${2:4}" -o "${2:4}" -s "ORDER BY chrom, strand, CASE strand WHEN '+' THEN txStart WHEN '-' THEN txEnd END"
        atdp "$@"
    inputBinding:
      position: 1
    doc: |
      Bash function to run samtools sort with all input parameters or skip it if trigger is false

  input_file:
    type:
      - File
    inputBinding:
      position: 2
      prefix: --in=
      separate: false
    doc: |
      Input indexed BAM file (optionally +BAI index file in secondaryFiles)

  annotation_filename:
    type:
      - File
    inputBinding:
      position: 3
      prefix: --a=
      separate: false
    doc: |
      Annotation file, tsv

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 4
      prefix: --out=
      separate: false
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename('_atdp.tsv');
            } else {
              return self;
            }
        }
    default: ""
    doc: |
      Base output file name, tsv

  log_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 5
      prefix: --log=
      separate: false
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename('_atdp.log');
            } else {
              return self;
            }
        }
    default: ""
    doc: |
      Log filename

  fragmentsize_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 6
      prefix: --f=
      separate: false
    doc: |
      Fragmentsize, int [150]

  avd_window_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 7
      prefix: --avd_window=
      separate: false
    doc: |
      Average tag density window, int [5000]

  avd_smooth_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 8
      prefix: --avd_smooth=
      separate: false
    doc: |
      Average smooth window (odd), int [0]

  ignore_chr:
    type:
      - "null"
      - string
    inputBinding:
      position: 9
      prefix: --sam_ignorechr=
      separate: false
    doc: |
      The chromosomes to be ignored, string

  double_chr:
    type:
      - "null"
      - string
    inputBinding:
      position: 10
      prefix: --sam_twicechr=
      separate: false
    doc: |
      The chromosomes to be doubled, string

  avd_heat_window_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 11
      prefix: --avd_heat_window=
      separate: false
    doc: |
      Average tag density window for heatmap

  mapped_reads:
    type:
      - "null"
      - int
    inputBinding:
      position: 12
      prefix: --m=
      separate: false
    doc: |
      Mapped reads number

  index_file:
    type:
      - "null"
      - File
    inputBinding:
      position: 13
      prefix: --index=
      separate: false
    doc: |
      Index file

outputs:
  log_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.log_filename == ""){
            return default_output_filename('_atdp.log');
          } else {
            return inputs.log_filename;
          }
        }

  result_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == ""){
            return default_output_filename('_atdp.tsv');
          } else {
            return inputs.output_filename;
          }
        }

baseCommand: [bash, '-c']