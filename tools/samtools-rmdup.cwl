#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: ./metadata/envvar-global.yml
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
            if (self == null || inputs.trigger == false){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Writes the output bam file to output_filename if set,
      otherwise generates output_filename on the base of bam_file

outputs:
  rmdup_output:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null || inputs.trigger == false){
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
          if (inputs.output_filename == null || inputs.trigger == false){
            return default_output_filename().split('.').slice(0,-1).join('.') + '.rmdup';
          } else {
            return inputs.output_filename.split('.').slice(0,-1).join('.') + '.rmdup';
          }
        }

baseCommand: [bash, '-c']

arguments:
  - valueFrom: |
      ${
        if (inputs.output_filename == null || inputs.trigger == false){
          return " > " + default_output_filename().split('.').slice(0,-1).join('.') + ".rmdup 2>&1";
        } else {
          return " > " + inputs.output_filename.split('.').slice(0,-1).join('.') + ".rmdup 2>&1";
        }
      }
    position: 100000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:name: "samtools-rmdup"
s:downloadUrl: https://github.com/SciDAP/workflows/blob/master/tools/samtools-rmdup.cwl
s:codeRepository: https://github.com/SciDAP/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool to remove duplicates from coordinate sorted BAM file set as input `bam_file`.
  If input `trigger` is set to `true` or isn't set at all (`true` is used by default), run `samtools rmdup`, return
  newly generated BAM file without duplicates and log as outputs `rmdup_output` and `rmdup_log`.
  If input `trigger` is set to `false`, return unchanged BAM and index files (if provided in secondaryFiles),
  previously staged into output directory, and log.

  Before execution `baseCommand` input BAM and index files (if provided in secondaryFiles) are staged into directory
  set as docker parameter `--workdir` (tool's output directory), using `InitialWorkDirRequirement`. Setting
  `writable: true` makes cwl-runner to make copies of the BAM and index files (if provided in secondaryFiles) and mount
  them to docker container with `rw` mode as part of `--workdir` (if set to false, the files staged into output
  directory will be mounted to docker container separately with `ro` mode). Because `bam_file` is copied into `--workdir`
  (where `samtools rmdup` should be run) we need to change its name - set `entryname`. This prevents `samtools rmdup`
  from reading and writing to the same file (which happens only when `output_filename` is not set and we generated
  output filename automatically on the base of the `bam_file` basename). Index file which is optionally set as
  `secondaryFiles` for `bam_file` shouln't be renamed, because `samtools rmdup` don't work with this file at all and is
  completely ignored when `trigger` is true.

  Trigger logic is implemented in bash script set by default in input `bash_script`. If first argment $0 (which is `trigger` input)
  is true, run `samtools rmdup` with the rest of the arguments. If $0 is not true, skip samtools running and rename
  staged into output directory BAM file. Current filename of staged BAM file and the filename which we are waiting to be the
  output of tool are always set as two last arguments of the script.

  Input `trigger` is Boolean, but returns String, because of `valueFrom` field. The `valueFrom` is used, because if `trigger`
  is false, cwl-runner doesn't append this argument at all to the the `baseCommand` - new feature of CWL v1.0.2. Alternatively,
  `prefix` field could be used, but it causes changing logic in bash script saved in `bash_script` input.

  `default_output_filename` function is used for generating output filename if input `output_filename` is not set or in
  case when `trigger` is false and we need to return BAM and index files (if provided in secondaryFiles) staged into
  output directory.

  Output `rmdup_output` returns `secondaryFiles` only in case when `trigger` was set to false (we need to rerun index).
    
s:about: |
  Usage:  samtools rmdup [-sS] <input.srt.bam> <output.bam>

  Option: -s    rmdup for SE reads
          -S    treat PE reads as SE in rmdup (force -s)
        --input-fmt-option OPT[=VAL]
                 Specify a single input file format option in the form
                 of OPTION or OPTION=VALUE
        --output-fmt FORMAT[,OPT[=VAL]]...
                 Specify output format (SAM, BAM, CRAM)
        --output-fmt-option OPT[=VAL]
                 Specify a single output file format option in the form
                 of OPTION or OPTION=VALUE
        --reference FILE
                 Reference sequence FASTA FILE [null]
