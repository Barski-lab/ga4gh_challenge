cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
    - var default_output_filename = function() {
          return inputs.unsorted_file.location.split('/').slice(-1)[0];
      };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2

inputs:
  unsorted_file:
    type: File
    inputBinding:
      position: 4

  key:
    type:
      type: array
      items: string
      inputBinding:
        prefix: "-k"
    inputBinding:
      position: 1
    doc: |
      -k, --key=POS1[,POS2]
      start a key at POS1, end it at POS2 (origin 1)

  output_filename:
    type:
    - "null"
    - string
    doc: |
      Name for generated output file

outputs:
  sorted_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == null){
            return default_output_filename();
          } else {
            return inputs.output_filename;
          }
        }
    doc: |
      Sorted file

stdout: |
  ${
    if (inputs.output_filename == null){
      return default_output_filename();
    } else {
      return inputs.output_filename;
    }
  }

baseCommand: ["sort"]


