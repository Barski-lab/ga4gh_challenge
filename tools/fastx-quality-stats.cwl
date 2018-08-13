cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.') + ".fastxstat"
        };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastx_toolkit:v0.0.14


inputs:

  input_file:
    type: File
    inputBinding:
      position: 10
      prefix: -i
    doc: |
      FASTA/Q input file. If FASTA file is given, only nucleotides distribution is calculated (there's no quality info)

  new_output_format:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 5
      prefix: '-N'
    doc: |
      New output format (with more information per nucleotide/cycle).
      cycle (previously called 'column') = cycle number
      max-count
      For each nucleotide in the cycle (ALL/A/C/G/T/N):
          count   = number of bases found in this column.
          min     = Lowest quality score value found in this column.
          max     = Highest quality score value found in this column.
          sum     = Sum of quality score values for this column.
          mean    = Mean quality score value for this column.
          Q1	= 1st quartile quality score.
          med	= Median quality score.
          Q3	= 3rd quartile quality score.
          IQR	= Inter-Quartile range (Q3-Q1).
          lW	= 'Left-Whisker' value (for boxplotting).
          rW	= 'Right-Whisker' value (for boxplotting).

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 11
      prefix: -o
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: ""
    doc: |
      Output file to store generated statistics. If not provided - return from default_output_filename function

outputs:
  statistics_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == ""){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }
    doc: Generated statistics file

baseCommand: [fastx_quality_stats]
