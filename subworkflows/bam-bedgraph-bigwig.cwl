cwlVersion: v1.0
class: Workflow

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
  bam_file:
    type: File
    label: "Input BAM file"
    doc: "Input BAM file, sorted by coordinates"

  chrom_length_file:
    type: File
    label: "Chromosome length file"
    doc: "Tab delimited chromosome length file: <chromName><TAB><chromSize>"

  scale:
    type: float?
    label: "Genome coverage scaling coefficient"
    doc: "Coefficient to scale the genome coverage by a constant factor"

  mapped_reads_number:
    type: int?
    label: "Mapped reads number"
    doc: |
      Parameter to calculate scale as 1000000/mapped_reads_number. Ignored by bedtools-genomecov.cwl in
      bam_to_bedgraph step if scale is provided

  pairchip:
    type: boolean?
    label: "Enable paired-end genome coverage calculation"
    doc: "Enable paired-end genome coverage calculation"

  fragment_size:
    type: int?
    label: "Fixed fragment size"
    doc: "Set fixed fragment size for genome coverage calculation"

  strand:
    type: string?
    label: "Enable strand specific genome coverage calculation"
    doc: "Calculate genome coverage of intervals from a specific strand"

  bigwig_filename:
    type: string?
    label: "bigWig output filename"
    doc: "Output filename for generated bigWig"

  bedgraph_filename:
    type: string?
    label: "bedGraph output filename"
    doc: "Output filename for generated bedGraph"

  split:
    type: boolean?
    label: "Split reads by 'N' and 'D'"
    doc: "Calculate genome coverage for each part of the splitted by 'N' and 'D' read"

  dutp:
    type: boolean?
    label: "Enable dUTP"
    doc: "Change strand af the mate read, so both reads come from the same strand"


outputs:
  bigwig_file:
    type: File
    outputSource: sorted_bedgraph_to_bigwig/bigwig_file
    label: "bigWig output file"
    doc: "bigWig output file"

  bedgraph_file:
    type: File
    outputSource: sort_bedgraph/sorted_file
    label: "bedGraph output file"
    doc: "bedGraph output file"

steps:
  bam_to_bedgraph:
    run: ../tools/bedtools-genomecov.cwl
    in:
      input_file: bam_file
      depth:
        default: "-bg"
      split:
        source: split
        default: true
      output_filename: bedgraph_filename
      pairchip: pairchip
      fragment_size: fragment_size
      scale: scale
      mapped_reads_number: mapped_reads_number
      strand: strand
      du: dutp
    out: [genome_coverage_file]

  sort_bedgraph:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: bam_to_bedgraph/genome_coverage_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  sorted_bedgraph_to_bigwig:
    run: ../tools/ucsc-bedgraphtobigwig.cwl
    in:
      bedgraph_file: sort_bedgraph/sorted_file
      chrom_length_file: chrom_length_file
      output_filename: bigwig_filename
    out: [bigwig_file]