{
    "cwlVersion": "v1.0", 
    "$schemas": [
        "http://schema.org/docs/schema_org_rdfa.html"
    ], 
    "$graph": [
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "class": "EnvVarRequirement", 
                    "envDef": [
                        {
                            "envName": "PATH", 
                            "envValue": "/usr/local/bin/:/usr/bin:/bin"
                        }
                    ], 
                    "id": "#envvar-global.yml", 
                    "name": "#envvar-global.yml"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "ShellCommandRequirement"
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/bamtools:v2.4.1", 
                    "dockerFile": "$import: ../dockerfiles/bamtools/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": [
                        "File", 
                        {
                            "type": "array", 
                            "items": "File", 
                            "inputBinding": {
                                "prefix": "-in"
                            }
                        }
                    ], 
                    "inputBinding": {
                        "position": 2, 
                        "valueFrom": "${\n  if ( Object.prototype.toString.call(inputs.input_files) === '[object Array]'){\n    return null;\n  } else {\n    return [\"-in\", inputs.input_files.path];\n  }\n}\n"
                    }, 
                    "doc": "the input BAM file[s]\nNOTE: If cwl fix a bug https://github.com/common-workflow-language/common-workflow-language/issues/330\nwe'll be able to use MultipleInputFeatureRequirement for single-item array and it will work\neven without additional File type and valueFrom field\n", 
                    "id": "#bamtools-stats.cwl/input_files"
                }
            ], 
            "outputs": [
                {
                    "type": "double", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "stats.log", 
                        "outputEval": "${\n  var s = self[0].contents.replace(/ /g,'').replace(/ *\\([^)]*\\) */g,'');\n  var duplicates = parseInt(s.substring ( s.indexOf(\"Duplicates\")+11, s.indexOf(\"\\t\", (s.indexOf(\"Duplicates\")))  ));\n  return duplicates;\n}\n"
                    }, 
                    "id": "#bamtools-stats.cwl/duplicates"
                }, 
                {
                    "type": "double", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "stats.log", 
                        "outputEval": "${\n  var s = self[0].contents.replace(/ /g,'').replace(/ *\\([^)]*\\) */g,'');\n  var failedQC = parseInt(s.substring ( s.indexOf(\"FailedQC\")+9, s.indexOf(\"\\t\", (s.indexOf(\"FailedQC\")))  ));\n  return failedQC;\n}\n"
                    }, 
                    "id": "#bamtools-stats.cwl/failedQC"
                }, 
                {
                    "type": "double", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "stats.log", 
                        "outputEval": "${\n  var s = self[0].contents.replace(/ /g,'').replace(/ *\\([^)]*\\) */g,'');\n  var forwardstrand = parseInt(s.substring ( s.indexOf(\"Forwardstrand\")+14, s.indexOf(\"\\t\", (s.indexOf(\"Forwardstrand\")))  ));\n  return forwardstrand;\n}\n"
                    }, 
                    "id": "#bamtools-stats.cwl/forwardstrand"
                }, 
                {
                    "type": "double", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "stats.log", 
                        "outputEval": "${\n  var s = self[0].contents.replace(/ /g,'').replace(/ *\\([^)]*\\) */g,'');\n  var mappedreads = parseInt(s.substring ( s.indexOf(\"Mappedreads\")+12, s.indexOf(\"\\t\", (s.indexOf(\"Mappedreads\")))  ));\n  return mappedreads;\n}\n"
                    }, 
                    "id": "#bamtools-stats.cwl/mappedreads"
                }, 
                {
                    "type": "double", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "stats.log", 
                        "outputEval": "${\n  var s = self[0].contents.replace(/ /g,'').replace(/ *\\([^)]*\\) */g,'');\n  var pairedendreads = parseInt(s.substring ( s.indexOf(\"Paired-endreads\")+16, s.indexOf(\"\\t\", (s.indexOf(\"Paired-endreads\")))  ));\n  return pairedendreads;\n}\n"
                    }, 
                    "id": "#bamtools-stats.cwl/pairedendreads"
                }, 
                {
                    "type": "double", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "stats.log", 
                        "outputEval": "${\n  var s = self[0].contents.replace(/ /g,'').replace(/ *\\([^)]*\\) */g,'');\n  var reversestrand = parseInt(s.substring ( s.indexOf(\"Reversestrand\")+14, s.indexOf(\"\\t\", (s.indexOf(\"Reversestrand\")))  ));\n  return reversestrand;\n}\n"
                    }, 
                    "id": "#bamtools-stats.cwl/reversestrand"
                }, 
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "stats.log"
                    }, 
                    "id": "#bamtools-stats.cwl/stats_log"
                }, 
                {
                    "type": "double", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "stats.log", 
                        "outputEval": "${\n  var s = self[0].contents.replace(/ /g,'').replace(/ *\\([^)]*\\) */g,'');\n  var totalReads = parseInt(s.substring ( s.indexOf(\"Totalreads\")+11, s.indexOf(\"\\t\", (s.indexOf(\"Totalreads\")))  ));\n  return totalReads;\n}\n"
                    }, 
                    "id": "#bamtools-stats.cwl/totalReads"
                }
            ], 
            "baseCommand": [
                "bamtools"
            ], 
            "arguments": [
                "stats", 
                {
                    "valueFrom": "$('> ' + 'stats.log')", 
                    "position": 1000, 
                    "shellQuote": false
                }
            ], 
            "s:name": "bamtools-stats", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "Tool is used to calculate general alignment statistics from the input BAM file\n", 
            "id": "#bamtools-stats.cwl", 
            "http://schema.org/mainEntity": {
                "class": "http://schema.org/SoftwareSourceCode", 
                "s:about": "Software suite for programmers and end users that facilitates research analysis and data management using BAM files. BamTools provides both the first C++ API publicly available for BAM file support as well as a command-line toolkit.\n", 
                "s:codeRepository": "https://github.com/pezmaster31/bamtools", 
                "s:targetProduct": {
                    "class": "http://schema.org/SoftwareApplication", 
                    "s:applicationCategory": "command line tool", 
                    "http://schema.org/softwareVersion": "2.4.1"
                }, 
                "s:publication": [
                    {
                        "class": "http://schema.org/ScholarlyArticle", 
                        "id": "#bamtools-metadata.yaml/10.1186/gb-2009-10-3-r25", 
                        "s:author": [
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Derek W. Barnett"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Erik K. Garrison"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Aaron R. Quinlan"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Michael P. Str\u00f6mberg"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Gabor T. Marth"
                            }
                        ], 
                        "datePublished": "14 April 2011", 
                        "http://schema.org/name": "BamTools: a C++ API and toolkit for analyzing and managing BAM files", 
                        "http://schema.org/url": "https://academic.oup.com/bioinformatics/article/27/12/1691/255399/BamTools-a-C-API-and-toolkit-for-analyzing-and"
                    }
                ], 
                "s:creator": [
                    {
                        "class": "http://schema.org/Organization", 
                        "s:member": [
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "Main author", 
                                "http://schema.org/name": "Derek Barnett"
                            }
                        ], 
                        "http://schema.org/name": "Pacific Biosciences"
                    }
                ], 
                "id": "#bamtools-metadata.yaml", 
                "name": "#bamtools-metadata.yaml", 
                "http://schema.org/name": "bowtie2", 
                "http://schema.org/url": "https://github.com/pezmaster31/bamtools/wiki", 
                "http://schema.org/license": [
                    "https://opensource.org/licenses/MIT"
                ], 
                "http://schema.org/programmingLanguage": "C++", 
                "http://schema.org/discussionUrl": [
                    "https://github.com/pezmaster31/bamtools/issues"
                ]
            }, 
            "http://schema.org/downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/bamtools-stats.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:michael.kotliar@cchmc.org", 
                                            "http://schema.org/name": "Michael Kotliar", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "id": "#0000-0002-6486-3898"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "Usage: bamtools stats [-in <filename> -in <filename> ... | -list <filelist>] [statsOptions]\nInput & Output:\n  -in <BAM filename>                the input BAM file [stdin]\n  -list <filename>                  the input BAM file list, one\n                                    line per file\n\nAdditional Stats:\n  -insert                           summarize insert size data\n\nHelp:\n  --help, -h                        shows this help text\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/bedtools2:v2.26.0", 
                    "dockerFile": "$import: ../dockerfiles/bedtools/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": {
                        "name": "#bedtools-genomecov.cwl/dept/JustDepts", 
                        "type": "enum", 
                        "symbols": [
                            "#bedtools-genomecov.cwl/dept/JustDepts/-bg", 
                            "#bedtools-genomecov.cwl/dept/JustDepts/-bga", 
                            "#bedtools-genomecov.cwl/dept/JustDepts/-d", 
                            "#bedtools-genomecov.cwl/dept/JustDepts/-dz"
                        ]
                    }, 
                    "inputBinding": {
                        "position": 5
                    }, 
                    "id": "#bedtools-genomecov.cwl/dept"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "Change strand af the mate read (so both reads from the same strand) useful for strand specific.\nWorks for BAM files only\n", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-du"
                    }, 
                    "id": "#bedtools-genomecov.cwl/du"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "fixed fragment size", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-fs"
                    }, 
                    "id": "#bedtools-genomecov.cwl/fragmentsize"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "Input genome file.", 
                    "inputBinding": {
                        "position": 11, 
                        "prefix": "-g"
                    }, 
                    "id": "#bedtools-genomecov.cwl/genomeFile"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#bedtools-genomecov.cwl/genomecoverageout"
                }, 
                {
                    "type": "File", 
                    "doc": "The input file can be in BAM format\n    (Note: BAM _must_ be sorted by position)\nor <bed/gff/vcf>\n", 
                    "inputBinding": {
                        "position": 10, 
                        "valueFrom": "${\n  var prefix = ((/.*\\.bam$/i).test(inputs.input.path))?'-ibam':'-i';\n  return [prefix,inputs.input.path];\n}\n"
                    }, 
                    "secondaryFiles": "${\n if ((/.*\\.bam$/i).test(self.location))\n    return {\"location\": self.location+\".bai\", \"class\": \"File\"};\n return [];\n}\n", 
                    "id": "#bedtools-genomecov.cwl/input"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "Calculate coverage of 3\" positions (instead of entire interval).\n", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-3"
                    }, 
                    "id": "#bedtools-genomecov.cwl/m3"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "Calculate coverage of 5\" positions (instead of entire interval).\n", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-5"
                    }, 
                    "id": "#bedtools-genomecov.cwl/m5"
                }, 
                {
                    "type": [
                        "null", 
                        "double"
                    ], 
                    "doc": "Optional parameter to calculate scale as 1000000/mappedreads\n", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-scale", 
                        "valueFrom": "${\n  if (inputs.scale){\n    return null;\n  } else {\n    return 1000000/inputs.mappedreads;\n  }\n}\n"
                    }, 
                    "id": "#bedtools-genomecov.cwl/mappedreads"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "Combine all positions with a depth >= max into\na single bin in the histogram. Irrelevant\nfor -d and -bedGraph\n- (INTEGER)\n", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-max"
                    }, 
                    "id": "#bedtools-genomecov.cwl/max"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "pair-end chip seq experiment", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-pc"
                    }, 
                    "id": "#bedtools-genomecov.cwl/pairchip"
                }, 
                {
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "doc": "Scale the coverage by a constant factor.\nEach coverage value is multiplied by this factor before being reported.\nUseful for normalizing coverage by, e.g., reads per million (RPM).\n- Default is 1.0; i.e., unscaled.\n- (FLOAT)\n", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-scale"
                    }, 
                    "id": "#bedtools-genomecov.cwl/scale"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "treat \"split\" BAM or BED12 entries as distinct BED intervals.\nwhen computing coverage.\nFor BAM files, this uses the CIGAR \"N\" and \"D\" operations\nto infer the blocks for computing coverage.\nFor BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds\nfields (i.e., columns 10,11,12).\n", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-split"
                    }, 
                    "id": "#bedtools-genomecov.cwl/split"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "doc": "Calculate coverage of intervals from a specific strand.\nWith BED files, requires at least 6 columns (strand is column 6).\n- (STRING): can be + or -\n", 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-strand"
                    }, 
                    "id": "#bedtools-genomecov.cwl/strand"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "doc": "The file containing the genome coverage", 
                    "outputBinding": {
                        "glob": "${\n  if (inputs.genomecoverageout == null){\n    return inputs.input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".bed\";\n  } else {\n    return inputs.genomecoverageout;\n  }\n}\n"
                    }, 
                    "id": "#bedtools-genomecov.cwl/genomecoverage"
                }
            ], 
            "stdout": "${\n  if (inputs.genomecoverageout == null){\n    return inputs.input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".bed\";\n  } else {\n    return inputs.genomecoverageout;\n  }\n}\n", 
            "baseCommand": [
                "bedtools", 
                "genomecov"
            ], 
            "s:name": "bedtools-genomecov", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "Tool is used to calculate statistics on the base of FASTQ file quality scores\n", 
            "id": "#bedtools-genomecov.cwl", 
            "http://schema.org/mainEntity": {
                "class": "http://schema.org/SoftwareSourceCode", 
                "s:about": "A software suite for the comparison, manipulation and annotation of genomic features in browser extensible data (BED) and general feature format (GFF) format. BEDTools also supports the comparison of sequence alignments in BAM format to both BED and GFF features. The tools are extremely efficient and allow the user to compare large datasets (e.g. next-generation sequencing data) with both public and custom genome annotation tracks. BEDTools can be combined with one another as well as with standard UNIX commands, thus facilitating routine genomics tasks as well as pipelines that can quickly answer intricate questions of large genomic datasets.\n", 
                "s:codeRepository": "https://github.com/arq5x/bedtools2", 
                "s:targetProduct": {
                    "class": "http://schema.org/SoftwareApplication", 
                    "s:applicationCategory": "command line tool", 
                    "http://schema.org/softwareVersion": "2.26.0"
                }, 
                "s:publication": [
                    {
                        "class": "http://schema.org/ScholarlyArticle", 
                        "id": "#bedtools-metadata.yaml/10.1093/bioinformatics/btq033", 
                        "s:author": [
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Aaron R. Quinlan"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Ira M. Hall"
                            }
                        ], 
                        "datePublished": "28 January 2010", 
                        "http://schema.org/name": "BEDTools: a flexible suite of utilities for comparing genomic features", 
                        "http://schema.org/url": "https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq033"
                    }
                ], 
                "s:creator": [
                    {
                        "class": "http://schema.org/Organization", 
                        "s:member": [
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "Main author", 
                                "http://schema.org/name": "Aaron Quinlan"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "Main author", 
                                "http://schema.org/name": "Neil Kindlon"
                            }
                        ], 
                        "http://schema.org/name": "University of Utah"
                    }
                ], 
                "id": "#bedtools-metadata.yaml", 
                "name": "#bedtools-metadata.yaml", 
                "http://schema.org/name": "bowtie2", 
                "http://schema.org/url": "http://bedtools.readthedocs.org", 
                "http://schema.org/license": [
                    "https://opensource.org/licenses/GPL-2.0"
                ], 
                "http://schema.org/programmingLanguage": "C++", 
                "http://schema.org/discussionUrl": [
                    "https://github.com/arq5x/bedtools2/issues", 
                    "https://groups.google.com/forum/#!forum/bedtools-discuss"
                ]
            }, 
            "http://schema.org/downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/bedtools-genomecov.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:Andrey.Kartashov@cchmc.org", 
                                            "http://schema.org/name": "Andrey Kartashov", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "id": "#0000-0001-9102-5681"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "Usage: bedtools genomecov [OPTIONS] -i <bed/gff/vcf> -g <genome>\nOptions:\n\t-ibam\t\tThe input file is in BAM format.\n\t\t\tNote: BAM _must_ be sorted by position\n\n\t-d\t\tReport the depth at each genome position (with one-based coordinates).\n\t\t\tDefault behavior is to report a histogram.\n\n\t-dz\t\tReport the depth at each genome position (with zero-based coordinates).\n\t\t\tReports only non-zero positions.\n\t\t\tDefault behavior is to report a histogram.\n\n\t-bg\t\tReport depth in BedGraph format. For details, see:\n\t\t\tgenome.ucsc.edu/goldenPath/help/bedgraph.html\n\n\t-bga\t\tReport depth in BedGraph format, as above (-bg).\n\t\t\tHowever with this option, regions with zero\n\t\t\tcoverage are also reported. This allows one to\n\t\t\tquickly extract all regions of a genome with 0\n\t\t\tcoverage by applying: \"grep -w 0$\" to the output.\n\n\t-split\t\tTreat \"split\" BAM or BED12 entries as distinct BED intervals.\n\t\t\twhen computing coverage.\n\t\t\tFor BAM files, this uses the CIGAR \"N\" and \"D\" operations\n\t\t\tto infer the blocks for computing coverage.\n\t\t\tFor BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds\n\t\t\tfields (i.e., columns 10,11,12).\n\n\t-strand\t\tCalculate coverage of intervals from a specific strand.\n\t\t\tWith BED files, requires at least 6 columns (strand is column 6).\n\t\t\t- (STRING): can be + or -\n\n\t-pc\t\tCalculate coverage of pair-end fragments.\n\t\t\tWorks for BAM files only\n\t-fs\t\tForce to use provided fragment size instead of read length\n\t\t\tWorks for BAM files only\n\t-du\t\tChange strand af the mate read (so both reads from the same strand) useful for strand specific\n\t\t\tWorks for BAM files only\n\t-5\t\tCalculate coverage of 5\" positions (instead of entire interval).\n\n\t-3\t\tCalculate coverage of 3\" positions (instead of entire interval).\n\n\t-max\t\tCombine all positions with a depth >= max into\n\t\t\ta single bin in the histogram. Irrelevant\n\t\t\tfor -d and -bedGraph\n\t\t\t- (INTEGER)\n\n\t-scale\t\tScale the coverage by a constant factor.\n\t\t\tEach coverage value is multiplied by this factor before being reported.\n\t\t\tUseful for normalizing coverage by, e.g., reads per million (RPM).\n\t\t\t- Default is 1.0; i.e., unscaled.\n\t\t\t- (FLOAT)\n\n\t-trackline\tAdds a UCSC/Genome-Browser track line definition in the first line of the output.\n\t\t\t- See here for more details about track line definition:\n\t\t\t      http://genome.ucsc.edu/goldenPath/help/bedgraph.html\n\t\t\t- NOTE: When adding a trackline definition, the output BedGraph can be easily\n\t\t\t      uploaded to the Genome Browser as a custom track,\n\t\t\t      BUT CAN NOT be converted into a BigWig file (w/o removing the first line).\n\n\t-trackopts\tWrites additional track line definition parameters in the first line.\n\t\t\t- Example:\n\t\t\t   -trackopts 'name=\"My Track\" visibility=2 color=255,30,30'\n\t\t\t   Note the use of single-quotes if you have spaces in your parameters.\n\t\t\t- (TEXT)\n\nNotes:\n\t(1) The genome file should tab delimited and structured as follows:\n\t <chromName><TAB><chromSize>\n\n\tFor example, Human (hg19):\n\tchr1\t249250621\n\tchr2\t243199373\n\t...\n\tchr18_gl000207_random\t4262\n\n\t(2) The input BED (-i) file must be grouped by chromosome.\n\t A simple \"sort -k 1,1 <BED> > <BED>.sorted\" will suffice.\n\n\t(3) The input BAM (-ibam) file must be sorted by position.\n\t A \"samtools sort <BAM>\" should suffice.\n\nTips:\n\tOne can use the UCSC Genome Browser's MySQL database to extract\n\tchromosome sizes. For example, H. sapiens:\n\n\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \\\n\t\"select chrom, size from hg19.chromInfo\" > hg19.genome\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "ShellCommandRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement", 
                    "expressionLib": [
                        "var default_output_filename = function() { if (Array.isArray(inputs.filelist) && inputs.filelist.length > 0){ return inputs.filelist[0].location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".sam\"; } else if (inputs.filelist != null){ return inputs.filelist.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".sam\"; } else if (Array.isArray(inputs.filelist_mates) && inputs.filelist_mates.length > 0){ return inputs.filelist_mates[0].location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".sam\"; } else if (inputs.filelist_mates != null){ return inputs.filelist_mates.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".sam\"; } else { return null; } };"
                    ]
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/bowtie:v1.2.0", 
                    "dockerFile": "$import: ../dockerfiles/bowtie/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-B"
                    }, 
                    "id": "#bowtie.cwl/B"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "reads and index are in colorspace\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-C"
                    }, 
                    "id": "#bowtie.cwl/C"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--minins <int>  minimum insert size for paired-end alignment (default: 0)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-I"
                    }, 
                    "id": "#bowtie.cwl/I"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>           like -m, but reports 1 random hit (MAPQ=0); requires --best\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-M"
                    }, 
                    "id": "#bowtie.cwl/M"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-Q"
                    }, 
                    "id": "#bowtie.cwl/Q"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--Q1"
                    }, 
                    "id": "#bowtie.cwl/Q1"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--maxins <int>  maximum insert size for paired-end alignment (default: 250)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-X"
                    }, 
                    "id": "#bowtie.cwl/X"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "--all           report all alignments per read (much slower than low -k)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-a"
                    }, 
                    "id": "#bowtie.cwl/a"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "<fname>       write aligned reads/pairs to file(s) <fname>\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--al"
                    }, 
                    "id": "#bowtie.cwl/al"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "hits guaranteed best stratum; ties broken by quality\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--best"
                    }, 
                    "id": "#bowtie.cwl/best"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "query sequences given on cmd line (as <mates>, <singles>)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-c"
                    }, 
                    "id": "#bowtie.cwl/c"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>   max megabytes of RAM for best-first search frames (def: 64)\nReporting:\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--chunkmbs"
                    }, 
                    "id": "#bowtie.cwl/chunkmbs"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--trim3 <int>   trim <int> bases from 3' (right) end of reads\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-3"
                    }, 
                    "id": "#bowtie.cwl/clip_3p_end"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--trim5 <int>   trim <int> bases from 5' (left) end of reads\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-5"
                    }, 
                    "id": "#bowtie.cwl/clip_5p_end"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "print original colorspace quals, not decoded quals\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--col-cqual"
                    }, 
                    "id": "#bowtie.cwl/col-cqual"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "print aligned colorspace seqs as colors, not decoded bases\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--col-cseq"
                    }, 
                    "id": "#bowtie.cwl/col-cseq"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "keep nucleotides at extreme ends of decoded alignment\nSAM:\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--col-keepends"
                    }, 
                    "id": "#bowtie.cwl/col-keepends"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-e"
                    }, 
                    "id": "#bowtie.cwl/e"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "query input files are (multi-)FASTA .fa/.mfa\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-f"
                    }, 
                    "id": "#bowtie.cwl/f"
                }, 
                {
                    "type": [
                        "null", 
                        "File", 
                        {
                            "type": "array", 
                            "items": "File"
                        }
                    ], 
                    "doc": "{-1 <m1> -2 <m2> | --12 <r> | <s>}\n<m1>    Comma-separated list of files containing upstream mates (or the\n      sequences themselves, if -c is set) paired with mates in <m2>\n<m2>    Comma-separated list of files containing downstream mates (or the\n      sequences themselves if -c is set) paired with mates in <m1>\n<r>     Comma-separated list of files containing Crossbow-style reads.  Can be\n      a mixture of paired and unpaired.  Specify \"-\"for stdin.\n<s>     Comma-separated list of files containing unpaired reads, or the\n      sequences themselves, if -c is set.  Specify \"-\"for stdin.\n", 
                    "inputBinding": {
                        "itemSeparator": ",", 
                        "position": 83
                    }, 
                    "id": "#bowtie.cwl/filelist"
                }, 
                {
                    "type": [
                        "null", 
                        {
                            "type": "array", 
                            "items": "File"
                        }
                    ], 
                    "inputBinding": {
                        "itemSeparator": ",", 
                        "position": 86, 
                        "prefix": "-12"
                    }, 
                    "doc": "Comma-separated list of files containing Crossbow-style reads.\nCan be a mixture of paired and unpaired.  Specify \"-\"for stdin.\n", 
                    "id": "#bowtie.cwl/filelist_crossbow"
                }, 
                {
                    "type": [
                        "null", 
                        "File", 
                        {
                            "type": "array", 
                            "items": "File"
                        }
                    ], 
                    "inputBinding": {
                        "itemSeparator": ",", 
                        "position": 85
                    }, 
                    "id": "#bowtie.cwl/filelist_mates"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 90, 
                        "valueFrom": "${\n    if (self == null){\n      return default_output_filename();\n    } else {\n      return self;\n    }\n}\n"
                    }, 
                    "default": null, 
                    "doc": "Generates default output filename on the base of filelist/filelist_mates files\n", 
                    "id": "#bowtie.cwl/filename"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--fr"
                    }, 
                    "id": "#bowtie.cwl/fr"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "write entire ref name (default: only up to 1st space)\nColorspace:\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--fullref"
                    }, 
                    "id": "#bowtie.cwl/fullref"
                }, 
                {
                    "type": "Directory", 
                    "doc": "Folder with indices files\n", 
                    "inputBinding": {
                        "position": 81, 
                        "valueFrom": "${\n    for (var i = 0; i < self.listing.length; i++) {\n        if (self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.ebwt' ||\n            self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.ebwtl'){\n          return self.listing[i].path.split('.').slice(0,-3).join('.');\n        }\n    }\n    return null;\n}\n"
                    }, 
                    "id": "#bowtie.cwl/indices_folder"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "qualities are given as space-separated integers (not ASCII)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--integer-quals"
                    }, 
                    "id": "#bowtie.cwl/integer-quals"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>           report up to <int> good alignments per read (default: 1)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-k"
                    }, 
                    "id": "#bowtie.cwl/k"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--seedlen <int> seed length for -n (default: 28)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-l"
                    }, 
                    "id": "#bowtie.cwl/l"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "force usage of a 'large' index, even if a small one is present\nAlignment:\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--large-index"
                    }, 
                    "id": "#bowtie.cwl/large-index"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>           suppress all alignments if > <int> exist (def: no limit)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-m"
                    }, 
                    "id": "#bowtie.cwl/m"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>       default mapping quality (MAPQ) to print for SAM alignments\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--mapq"
                    }, 
                    "id": "#bowtie.cwl/mapq"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "<fname>      write reads/pairs over -m limit to file(s) <fname>\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--max"
                    }, 
                    "id": "#bowtie.cwl/max"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--maxbts"
                    }, 
                    "id": "#bowtie.cwl/maxbts"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "use memory-mapped I/O for index; many 'bowtie's can share\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--mm"
                    }, 
                    "id": "#bowtie.cwl/mm"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-n"
                    }, 
                    "id": "#bowtie.cwl/n"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "--norc      do not align to forward/reverse-complement reference strand\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--nofw"
                    }, 
                    "id": "#bowtie.cwl/nofw"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "disable Maq-like quality rounding for -n (nearest 10 <= 30)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--nomaqround"
                    }, 
                    "id": "#bowtie.cwl/nomaqround"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--offrate <int> override offrate of index; must be >= index's offrate\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-o"
                    }, 
                    "id": "#bowtie.cwl/o"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>  max # attempts to find mate for anchor hit (default: 100)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--pairtries"
                    }, 
                    "id": "#bowtie.cwl/pairtries"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "input quals are Phred+33 (default)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--phred33-quals"
                    }, 
                    "id": "#bowtie.cwl/phred33-quals"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "input quals are Phred+64 (same as --solexa1.3-quals)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--phred64-quals"
                    }, 
                    "id": "#bowtie.cwl/phred64-quals"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "query input files are FASTQ .fq/.fastq (default)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-q"
                    }, 
                    "id": "#bowtie.cwl/q"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "print nothing but the alignments\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--quiet"
                    }, 
                    "id": "#bowtie.cwl/quiet"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "query input files are raw one-sequence-per-line\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-r"
                    }, 
                    "id": "#bowtie.cwl/r"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "# of reads to read from input file at once (default: 16)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--reads-per-batch"
                    }, 
                    "id": "#bowtie.cwl/reads_per_batch"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "refer to ref. seqs by 0-based index rather than name\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--refidx"
                    }, 
                    "id": "#bowtie.cwl/refidx"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "write alignments to files refXXXXX.map, 1 map per reference\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--refout"
                    }, 
                    "id": "#bowtie.cwl/refout"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--skip <int>    skip the first <int> reads/pairs in the input\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-s"
                    }, 
                    "id": "#bowtie.cwl/s"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "--sam           write hits in SAM format\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-S"
                    }, 
                    "id": "#bowtie.cwl/sam"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "doc": "<text>    add <text> (usually \"lab=value\") to @RG line of SAM header\nPerformance:\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--sam-RG"
                    }, 
                    "id": "#bowtie.cwl/sam-RG"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "supppress header lines (starting with @) for SAM output\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--sam-nohead"
                    }, 
                    "id": "#bowtie.cwl/sam-nohead"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "supppress @SQ header lines for SAM output\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--sam-nosq"
                    }, 
                    "id": "#bowtie.cwl/sam-nosq"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>       seed for random number generator\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--seed"
                    }, 
                    "id": "#bowtie.cwl/seed"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "use shared mem for index; many 'bowtie's can share\nOther:\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--shmem"
                    }, 
                    "id": "#bowtie.cwl/shmem"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "<dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--snpfrac"
                    }, 
                    "id": "#bowtie.cwl/snpfrac"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>   Phred penalty for SNP when decoding colorspace (def: 30)\nor\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--snpphred"
                    }, 
                    "id": "#bowtie.cwl/snpphred"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "input quals are from GA Pipeline ver. < 1.3\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--solexa-quals"
                    }, 
                    "id": "#bowtie.cwl/solexa-quals"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "input quals are from GA Pipeline ver. >= 1.3\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--solexa1.3-quals"
                    }, 
                    "id": "#bowtie.cwl/solexa1.3-quals"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "hits in sub-optimal strata aren't reported (requires --best)\nOutput:\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--strata"
                    }, 
                    "id": "#bowtie.cwl/strata"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "<cols>  suppresses given columns (comma-delim'ed) in default output\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--suppress"
                    }, 
                    "id": "#bowtie.cwl/suppress"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "--time          print wall-clock time taken by search phases\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-t"
                    }, 
                    "id": "#bowtie.cwl/t"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--threads <int> number of alignment threads to launch (default: 1)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-p"
                    }, 
                    "id": "#bowtie.cwl/threads"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-u"
                    }, 
                    "id": "#bowtie.cwl/u"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "<fname>       write unaligned reads/pairs to file(s) <fname>\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--un"
                    }, 
                    "id": "#bowtie.cwl/un"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "<int>           report end-to-end hits w/ <=v mismatches; ignore qualities\nor\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-v"
                    }, 
                    "id": "#bowtie.cwl/v"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "verbose output (for debugging)\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "--verbose"
                    }, 
                    "id": "#bowtie.cwl/verbose"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "--tryhard       try hard to find valid alignments, at the expense of speed\n", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-y"
                    }, 
                    "id": "#bowtie.cwl/y"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "${\n   if (inputs.filename == null){\n     return default_output_filename();\n   } else {\n     return inputs.filename;\n   }\n}\n"
                    }, 
                    "id": "#bowtie.cwl/output"
                }, 
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "${\n   if (inputs.filename == null){\n     return default_output_filename() + \".log\";\n   } else {\n     return inputs.filename + \".log\";\n   }\n}\n"
                    }, 
                    "id": "#bowtie.cwl/output_bowtie_log"
                }
            ], 
            "baseCommand": [
                "bowtie"
            ], 
            "arguments": [
                {
                    "valueFrom": "${\n  if (inputs.filelist && inputs.filelist_mates){\n    return \"-1\";\n  }\n  return null;\n}\n", 
                    "position": 82
                }, 
                {
                    "valueFrom": "${\n  if (inputs.filelist && inputs.filelist_mates){\n    return \"-2\";\n  }\n  return null;\n}\n", 
                    "position": 84
                }, 
                {
                    "valueFrom": "${\n  if (inputs.filename == null){\n    return ' 2> ' + default_output_filename() + '.log';\n  } else {\n    return ' 2> ' + inputs.filename + '.log';\n  }\n}\n", 
                    "position": 100000, 
                    "shellQuote": false
                }
            ], 
            "s:name": "bowtie", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "Tool is used to run bowtie aligner to align input FASTQ file(s) to reference genome\n", 
            "id": "#bowtie.cwl", 
            "http://schema.org/mainEntity": {
                "class": "http://schema.org/SoftwareSourceCode", 
                "s:about": "Bowtie is an ultrafast, memory-efficient short read aligner. It aligns short DNA sequences (reads) to the human genome at a rate of over 25 million 35-bp reads per hour. Bowtie indexes the genome with a Burrows-Wheeler index to keep its memory footprint small: typically about 2.2 GB for the human genome (2.9 GB for paired-end).\n", 
                "s:codeRepository": "https://github.com/BenLangmead/bowtie", 
                "s:targetProduct": {
                    "class": "http://schema.org/SoftwareApplication", 
                    "s:applicationCategory": "command line tool", 
                    "http://schema.org/softwareVersion": "1.2.0"
                }, 
                "s:publication": [
                    {
                        "class": "http://schema.org/ScholarlyArticle", 
                        "id": "#bowtie-metadata.yaml/10.1186/gb-2009-10-3-r25", 
                        "s:author": [
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Ben Langmead"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Steven L Salzberg"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Cole Trapnell"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Mihai Pop"
                            }
                        ], 
                        "datePublished": "04 March 2009", 
                        "http://schema.org/name": "Ultrafast and memory-efficient alignment of short DNA sequences to the human genome", 
                        "http://schema.org/url": "https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r25"
                    }
                ], 
                "s:creator": [
                    {
                        "class": "http://schema.org/Organization", 
                        "s:member": [
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "Main author", 
                                "http://schema.org/name": "Ben Langmead"
                            }
                        ], 
                        "http://schema.org/name": "Johns Hopkins University"
                    }
                ], 
                "id": "#bowtie-metadata.yaml", 
                "name": "#bowtie-metadata.yaml", 
                "http://schema.org/name": "bowtie2", 
                "http://schema.org/url": "http://bowtie-bio.sourceforge.net/index.shtml", 
                "http://schema.org/license": [
                    "https://opensource.org/licenses/GPL-3.0"
                ], 
                "http://schema.org/programmingLanguage": "C++", 
                "http://schema.org/discussionUrl": [
                    "https://github.com/BenLangmead/bowtie/issues"
                ]
            }, 
            "http://schema.org/downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/bowtie.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:Andrey.Kartashov@cchmc.org", 
                                            "http://schema.org/name": "Andrey Kartashov", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0001-9102-5681"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "Usage: bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]\n\n  <m1>    Comma-separated list of files containing upstream mates (or the\n          sequences themselves, if -c is set) paired with mates in <m2>\n  <m2>    Comma-separated list of files containing downstream mates (or the\n          sequences themselves if -c is set) paired with mates in <m1>\n  <r>     Comma-separated list of files containing Crossbow-style reads.  Can be\n          a mixture of paired and unpaired.  Specify \"-\" for stdin.\n  <s>     Comma-separated list of files containing unpaired reads, or the\n          sequences themselves, if -c is set.  Specify \"-\" for stdin.\n  <hit>   File to write hits to (default: stdout)\nInput:\n  -q                 query input files are FASTQ .fq/.fastq (default)\n  -f                 query input files are (multi-)FASTA .fa/.mfa\n  -r                 query input files are raw one-sequence-per-line\n  -c                 query sequences given on cmd line (as <mates>, <singles>)\n  -C                 reads and index are in colorspace\n  -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C\n  --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively\n  -s/--skip <int>    skip the first <int> reads/pairs in the input\n  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)\n  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads\n  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads\n  --phred33-quals    input quals are Phred+33 (default)\n  --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)\n  --solexa-quals     input quals are from GA Pipeline ver. < 1.3\n  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3\n  --integer-quals    qualities are given as space-separated integers (not ASCII)\n  --large-index      force usage of a 'large' index, even if a small one is present\nAlignment:\n  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities\n    or\n  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)\n  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)\n  -l/--seedlen <int> seed length for -n (default: 28)\n  --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)\n  -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)\n  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)\n  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)\n  --nofw/--norc      do not align to forward/reverse-complement reference strand\n  --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)\n  --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)\n  -y/--tryhard       try hard to find valid alignments, at the expense of speed\n  --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)\n --reads-per-batch   # of reads to read from input file at once (default: 16)\nReporting:\n  -k <int>           report up to <int> good alignments per read (default: 1)\n  -a/--all           report all alignments per read (much slower than low -k)\n  -m <int>           suppress all alignments if > <int> exist (def: no limit)\n  -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best\n  --best             hits guaranteed best stratum; ties broken by quality\n  --strata           hits in sub-optimal strata aren't reported (requires --best)\nOutput:\n  -t/--time          print wall-clock time taken by search phases\n  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)\n  --quiet            print nothing but the alignments\n  --refout           write alignments to files refXXXXX.map, 1 map per reference\n  --refidx           refer to ref. seqs by 0-based index rather than name\n  --al <fname>       write aligned reads/pairs to file(s) <fname>\n  --un <fname>       write unaligned reads/pairs to file(s) <fname>\n  --max <fname>      write reads/pairs over -m limit to file(s) <fname>\n  --suppress <cols>  suppresses given columns (comma-delim'ed) in default output\n  --fullref          write entire ref name (default: only up to 1st space)\nColorspace:\n  --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)\n     or\n  --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred\n  --col-cseq         print aligned colorspace seqs as colors, not decoded bases\n  --col-cqual        print original colorspace quals, not decoded quals\n  --col-keepends     keep nucleotides at extreme ends of decoded alignment\nSAM:\n  -S/--sam           write hits in SAM format\n  --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments\n  --sam-nohead       supppress header lines (starting with @) for SAM output\n  --sam-nosq         supppress @SQ header lines for SAM output\n  --sam-RG <text>    add <text> (usually \"lab=value\") to @RG line of SAM header\nPerformance:\n  -o/--offrate <int> override offrate of index; must be >= index's offrate\n  -p/--threads <int> number of alignment threads to launch (default: 1)\n  --mm               use memory-mapped I/O for index; many 'bowtie's can share\n  --shmem            use shared mem for index; many 'bowtie's can share\nOther:\n  --seed <int>       seed for random number generator\n  --verbose          verbose output (for debugging)\n  --version          print version information and quit\n  -h/--help          print this usage message\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/fastx_toolkit:v0.0.14", 
                    "dockerFile": "$import: ../dockerfiles/fastx/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 10, 
                        "prefix": "-i"
                    }, 
                    "doc": "FASTA/Q input file. If FASTA file is given, only nucleotides distribution is calculated (there's no quality info)\n", 
                    "id": "#fastx-quality-stats.cwl/input_file"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 5, 
                        "prefix": "-N"
                    }, 
                    "doc": "New output format (with more information per nucleotide/cycle).\ncycle (previously called 'column') = cycle number\nmax-count\nFor each nucleotide in the cycle (ALL/A/C/G/T/N):\n    count   = number of bases found in this column.\n    min     = Lowest quality score value found in this column.\n    max     = Highest quality score value found in this column.\n    sum     = Sum of quality score values for this column.\n    mean    = Mean quality score value for this column.\n    Q1\t= 1st quartile quality score.\n    med\t= Median quality score.\n    Q3\t= 3rd quartile quality score.\n    IQR\t= Inter-Quartile range (Q3-Q1).\n    lW\t= 'Left-Whisker' value (for boxplotting).\n    rW\t= 'Right-Whisker' value (for boxplotting).\n", 
                    "id": "#fastx-quality-stats.cwl/new_output_format"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "*.fastxstat"
                    }, 
                    "doc": "Statistics file", 
                    "id": "#fastx-quality-stats.cwl/statistics"
                }
            ], 
            "baseCommand": [
                "fastx_quality_stats"
            ], 
            "arguments": [
                {
                    "valueFrom": "$(inputs.input_file.basename + \".fastxstat\")", 
                    "position": 11, 
                    "prefix": "-o"
                }
            ], 
            "s:name": "fastx-quality-stats", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "Tool is used to calculate statistics on the base of FASTQ file quality scores\n", 
            "id": "#fastx-quality-stats.cwl", 
            "http://schema.org/mainEntity": {
                "class": "http://schema.org/SoftwareSourceCode", 
                "s:about": "The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing. Available Tools:\n  FASTQ-to-FASTA converter (Convert FASTQ files to FASTA files)\n  FASTQ Information (Chart Quality Statistics and Nucleotide Distribution)\n  FASTQ/A Collapser (Collapsing identical sequences in a FASTQ/A file into a single sequence)\n  FASTQ/A Trimmer (Shortening reads in a FASTQ or FASTQ files, removing barcodes or noise)\n  FASTQ/A Renamer (Renames the sequence identifiers in FASTQ/A file)\n  FASTQ/A Clipper (Removing sequencing adapters/linkers)\n  FASTQ/A Reverse-Complement (Producing the Reverse-complement of each sequence in a FASTQ/FASTA file)\n  FASTQ/A Barcode splitter (Splitting a FASTQ/FASTA files containning multiple samples)\n  FASTA Formatter (changes the width of sequences line in a FASTA file)\n  FASTA Nucleotide Changer (Convets FASTA sequences from/to RNA/DNA)\n  FASTQ Quality Filter (Filters sequences based on quality)\n  FASTQ Quality Trimmer (Trims (cuts) sequences based on quality)\n  FASTQ Masker (Masks nucleotides with 'N' (or other character) based on quality)\n", 
                "s:codeRepository": "https://github.com/agordon/fastx_toolkit", 
                "s:targetProduct": {
                    "class": "http://schema.org/SoftwareApplication", 
                    "s:applicationCategory": "command line tool", 
                    "http://schema.org/softwareVersion": "0.0.14"
                }, 
                "s:discussionUrl": [
                    "https://github.com/agordon/fastx_toolkit/issues"
                ], 
                "id": "#fastx-toolkit-metadata.yaml", 
                "name": "#fastx-toolkit-metadata.yaml", 
                "http://schema.org/name": "fastx-toolkit", 
                "http://schema.org/url": "http://hannonlab.cshl.edu/fastx_toolkit/index.html", 
                "http://schema.org/license": [
                    "https://opensource.org/licenses/GPL-3.0"
                ], 
                "http://schema.org/programmingLanguage": "C", 
                "http://schema.org/creator": [
                    {
                        "class": "http://schema.org/Organization", 
                        "s:member": [
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "Main author", 
                                "http://schema.org/name": "Assaf Gordon"
                            }
                        ], 
                        "http://schema.org/name": "The Hannon Lab"
                    }
                ]
            }, 
            "http://schema.org/downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/fastx-quality-stats.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:michael.kotliar@cchmc.org", 
                                            "http://schema.org/name": "Michael Kotliar", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0002-6486-3898"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "usage: fastx_quality_stats [-h] [-N] [-i INFILE] [-o OUTFILE] Part of FASTX Toolkit 0.0.14 by A. Gordon (assafgordon@gmail.com)\n   [-h] = This helpful help screen.\n   [-i INFILE]  = FASTQ input file. default is STDIN.\n   [-o OUTFILE] = TEXT output file. default is STDOUT.\n   [-N]         = New output format (with more information per nucleotide/cycle).\nThe *OLD* output TEXT file will have the following fields (one row per column):\n\tcolumn\t= column number (1 to 36 for a 36-cycles read solexa file)\n\tcount   = number of bases found in this column.\n\tmin     = Lowest quality score value found in this column.\n\tmax     = Highest quality score value found in this column.\n\tsum     = Sum of quality score values for this column.\n\tmean    = Mean quality score value for this column.\n\tQ1\t= 1st quartile quality score.\n\tmed\t= Median quality score.\n\tQ3\t= 3rd quartile quality score.\n\tIQR\t= Inter-Quartile range (Q3-Q1).\n\tlW\t= 'Left-Whisker' value (for boxplotting).\n\trW\t= 'Right-Whisker' value (for boxplotting).\n\tA_Count\t= Count of 'A' nucleotides found in this column.\n\tC_Count\t= Count of 'C' nucleotides found in this column.\n\tG_Count\t= Count of 'G' nucleotides found in this column.\n\tT_Count\t= Count of 'T' nucleotides found in this column.\n\tN_Count = Count of 'N' nucleotides found in this column.\n\tmax-count = max. number of bases (in all cycles)\nThe *NEW* output format:\n\tcycle (previously called 'column') = cycle number\n\tmax-count\n  For each nucleotide in the cycle (ALL/A/C/G/T/N):\n\t\tcount   = number of bases found in this column.\n\t\tmin     = Lowest quality score value found in this column.\n\t\tmax     = Highest quality score value found in this column.\n\t\tsum     = Sum of quality score values for this column.\n\t\tmean    = Mean quality score value for this column.\n\t\tQ1\t= 1st quartile quality score.\n\t\tmed\t= Median quality score.\n\t\tQ3\t= 3rd quartile quality score.\n\t\tIQR\t= Inter-Quartile range (Q3-Q1).\n\t\tlW\t= 'Left-Whisker' value (for boxplotting).\n\t\trW\t= 'Right-Whisker' value (for boxplotting).\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "InlineJavascriptRequirement", 
                    "expressionLib": [
                        "var default_output_filename = function() { return inputs.input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".sorted.bed\"; };"
                    ]
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "ubuntu:15.10"
                }
            ], 
            "inputs": [
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 4
                    }, 
                    "id": "#linux-sort.cwl/input"
                }, 
                {
                    "type": {
                        "type": "array", 
                        "items": "string", 
                        "inputBinding": {
                            "prefix": "-k"
                        }
                    }, 
                    "inputBinding": {
                        "position": 1
                    }, 
                    "doc": "-k, --key=POS1[,POS2]\nstart a key at POS1, end it at POS2 (origin 1)\n", 
                    "id": "#linux-sort.cwl/key"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "doc": "The sorted file", 
                    "outputBinding": {
                        "glob": "$(default_output_filename())"
                    }, 
                    "id": "#linux-sort.cwl/sorted"
                }
            ], 
            "stdout": "$(default_output_filename())", 
            "baseCommand": [
                "sort"
            ], 
            "s:downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/linux-sort.cwl", 
            "s:license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "s:creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:Andrey.Kartashov@cchmc.org", 
                                            "http://schema.org/name": "Andrey Kartashov", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0001-9102-5681"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "doc": "Tool is used to run sort command with input data\n", 
            "id": "#linux-sort.cwl", 
            "http://schema.org/name": "linux-sort", 
            "http://schema.org/codeRepository": "https://github.com/SciDAP/workflows", 
            "http://schema.org/isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "ShellCommandRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement", 
                    "expressionLib": [
                        "var default_name = function(input_staged, sufix) { input_staged = input_staged || false; sufix = sufix || \"\"; if (inputs.trigger == false && input_staged){ return input_staged.basename; } else { if (Object.prototype.toString.call(inputs.treatment) === '[object Array]'){ return inputs.treatment[0].location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+sufix; } else { return inputs.treatment.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+sufix; } } }"
                    ]
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": "${\n  var listing = []\n  if (inputs.peak_xls_file_staged){\n    listing.push(inputs.peak_xls_file_staged);\n  }\n  if (inputs.narrow_peak_file_staged){\n    listing.push(inputs.narrow_peak_file_staged);\n  }\n  if (inputs.broad_peak_file_staged){\n    listing.push(inputs.broad_peak_file_staged);\n  }\n  if (inputs.gapped_peak_file_staged){\n    listing.push(inputs.gapped_peak_file_staged);\n  }\n  if (inputs.peak_summits_file_staged){\n    listing.push(inputs.peak_summits_file_staged);\n  }\n  if (inputs.moder_r_file_staged){\n    listing.push(inputs.moder_r_file_staged);\n  }\n  if (inputs.treat_pileup_bdg_file_staged){\n    listing.push(inputs.treat_pileup_bdg_file_staged);\n  }\n  if (inputs.control_lambda_bdg_file_staged){\n    listing.push(inputs.control_lambda_bdg_file_staged);\n  }\n  return listing;\n}\n"
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/macs2:v2.1.1", 
                    "dockerFile": "$import: ../dockerfiles/macs2/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 17, 
                        "prefix": "--bdg"
                    }, 
                    "doc": "If this flag is on, MACS will store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files.\nThe bedGraph files will be stored in current directory named NAME+\u2019_treat_pileup.bdg\u2019 for treatment data, NAME+\u2019_control_lambda.bdg\u2019\nfor local lambda values from control, NAME+\u2019_treat_pvalue.bdg\u2019 for Poisson pvalue scores (in -log10(pvalue) form),\nand NAME+\u2019_treat_qvalue.bdg\u2019 for q-value scores from\nBenjamini\u2013Hochberg\u2013Yekutieli procedure <http://en.wikipedia.org/wiki/False_discovery_rate#Dependent_tests>\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/bdg"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 36, 
                        "prefix": "--broad"
                    }, 
                    "doc": "If set, MACS will try to call broad peaks by linking nearby highly enriched\nregions. The linking region is controlled by another cutoff through --linking-cutoff.\nThe maximum linking region length is 4 times of d from MACS.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/broad"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 37, 
                        "prefix": "--broad-cutoff"
                    }, 
                    "doc": "Cutoff for broad region. This option is not available\nunless --broad is set. If -p is set, this is a pvalue\ncutoff, otherwise, it's a qvalue cutoff.\nDEFAULT: 0.1\n", 
                    "id": "#macs2-callpeak.cwl/broad_cutoff"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "For staging in a case of trigger is set to false", 
                    "id": "#macs2-callpeak.cwl/broad_peak_file_staged"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 16, 
                        "prefix": "--buffer-size"
                    }, 
                    "doc": "Buffer size for incrementally increasing internal\narray size to store reads alignment information. In\nmost cases, you don't have to change this parameter.\nHowever, if there are large number of\nchromosomes/contigs/scaffolds in your alignment, it's\nrecommended to specify a smaller buffer size in order\nto decrease memory usage (but it will take longer time\nto read alignment files). Minimum memory requested for\nreading an alignment file is about # of CHROMOSOME *\nBUFFER_SIZE * 2 Bytes.\nDEFAULT: 100000\n", 
                    "id": "#macs2-callpeak.cwl/buffer_size"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 21, 
                        "prefix": "--bw"
                    }, 
                    "doc": "Band width for picking regions to compute  fragment  size.  This\nvalue is only used while building the shifting model\nDEFAULT: 300\n", 
                    "id": "#macs2-callpeak.cwl/bw"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 39, 
                        "prefix": "--call-summits"
                    }, 
                    "doc": "If set, MACS will use a more sophisticated signal processing approach to\nfind subpeak summits in each enriched peak region.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/call_summits"
                }, 
                {
                    "type": [
                        "null", 
                        "File", 
                        {
                            "type": "array", 
                            "items": "File"
                        }
                    ], 
                    "inputBinding": {
                        "position": 12, 
                        "prefix": "-c"
                    }, 
                    "doc": "The control or mock data file. Please follow the same direction as for -t/\u2013treatment.\n", 
                    "id": "#macs2-callpeak.cwl/control"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "For staging in a case of trigger is set to false", 
                    "id": "#macs2-callpeak.cwl/control_lambda_bdg_file_staged"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 38, 
                        "prefix": "--cutoff-analysis"
                    }, 
                    "doc": "While set, MACS2 will analyze number or total length of peaks that can be\ncalled by different p-value cutoff then output a summary table to help user\ndecide a better cutoff. The table will be saved in NAME_cutoff_analysis.txt\nfile. Note, minlen and maxgap may affect the results. WARNING: May take ~30\nfolds longer time to finish.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/cutoff_analysis"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 31, 
                        "prefix": "--down-sample"
                    }, 
                    "doc": "When set, random sampling method will scale down the bigger sample. By default,\nMACS uses linear scaling. Warning: This option will make your result unstable\nand irreproducible since each time, random reads would be selected. Consider\nto use ''randsample'' script instead. <not implmented>If used together with\n--SPMR, 1 million unique reads will be randomly picked.</not implemented> Caution:\ndue to the implementation, the final number of selected reads may not be as\nyou expected!\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/down_sample"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 26, 
                        "prefix": "--extsize"
                    }, 
                    "doc": "The arbitrary extension size in bp. When nomodel is  true, MACS will use\nthis value as fragment size to  extend each read towards 3'' end, then pile\nthem up.  It''s exactly twice the number of obsolete SHIFTSIZE.  In previous\nlanguage, each read is moved 5''->3''  direction to middle of fragment by 1/2\nd, then  extended to both direction with 1/2 d. This is  equivalent to say each\nread is extended towards 5''->3''  into a d size fragment.EXTSIZE\nand  SHIFT can be combined when necessary. Check SHIFT  option.\nDEFAULT: 200\n", 
                    "id": "#macs2-callpeak.cwl/extsize"
                }, 
                {
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "inputBinding": {
                        "position": 40, 
                        "prefix": "--fe-cutoff"
                    }, 
                    "doc": "When set, the value will be used to filter out peaks\nwith low fold-enrichment. Note, MACS2 use 1.0 as\npseudocount while calculating fold-enrichment.\nDEFAULT: 1.0\n", 
                    "id": "#macs2-callpeak.cwl/fe-cutoff"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 23, 
                        "prefix": "--fix-bimodal"
                    }, 
                    "doc": "Whether turn on the auto pair model process. If set, when MACS failed to\nbuild paired model, it will use the nomodel settings, the --exsize parameter\nto extend each tags towards 3'' direction. Not to use this automate fixation\nis a default behavior now.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/fix_bimodal"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 13, 
                        "prefix": "-f"
                    }, 
                    "doc": "{AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE}, --format\n{AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE} Format of tag file,\n\"AUTO\", \"BED\" or \"ELAND\" or \"ELANDMULTI\" or \"ELANDEXPORT\" or \"SAM\" or \"BAM\"\nor \"BOWTIE\" or \"BAMPE\". The default AUTO option will let MACS decide which format\nthe file is. Note that MACS can''t detect \"BAMPE\" or \"BEDPE\" format with \"AUTO\",\nand you have to implicitly specify the format for \"BAMPE\" and \"BEDPE\".\nDEFAULT: AUTO\n", 
                    "id": "#macs2-callpeak.cwl/format"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "For staging in a case of trigger is set to false", 
                    "id": "#macs2-callpeak.cwl/gapped_peak_file_staged"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 14, 
                        "prefix": "-g"
                    }, 
                    "doc": "It\u2019s the mappable genome size or effective genome size which is defined as the genome size which can be sequenced.\nBecause of the repetitive features on the chromsomes, the actual mappable genome size will be smaller than the\noriginal size, about 90% or 70% of the genome size. The default hs \u2013 2.7e9 is recommended for UCSC human hg18\nassembly. Here are all precompiled parameters for effective genome size:\n  hs:\t2.7e9\n  mm:\t1.87e9\n  ce:\t9e7\n  dm:\t1.2e8\nDEFAULT: hs\n", 
                    "id": "#macs2-callpeak.cwl/genome_size"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 15, 
                        "prefix": "--keep-dup"
                    }, 
                    "doc": "It controls the MACS behavior towards duplicate tags\nat the exact same location -- the same coordination\nand the same strand. The 'auto' option makes MACS\ncalculate the maximum tags at the exact same location\nbased on binomal distribution using 1e-5 as pvalue\ncutoff; and the 'all' option keeps every tags. If an\ninteger is given, at most this number of tags will be\nkept at the same location. The default is to keep one\ntag at the same location.\nDEFAULT: 1\n", 
                    "id": "#macs2-callpeak.cwl/keep_dup"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 35, 
                        "prefix": "--llocal"
                    }, 
                    "doc": "The large nearby region in basepairs to calculate\ndynamic lambda. This is used to capture the surround\nbias. If you set this to 0, MACS will skip llocal\nlambda calculation. *Note* that MACS will always\nperform a d-size local lambda calculation. The final\nlocal bias should be the maximum of the lambda value\nfrom d, slocal, and llocal size windows.\nDEFAULT: 10000.\n", 
                    "id": "#macs2-callpeak.cwl/llocal"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 22, 
                        "prefix": "-m", 
                        "valueFrom": "${\n  return self.replace(/\\s+/g, ' ').split(' ');\n}\n"
                    }, 
                    "doc": "Select the regions within MFOLD range of high-\nconfidence enrichment ratio against background to\nbuild model. Fold-enrichment in regions must be lower\nthan upper limit, and higher than the lower limit. Use\nas \"-m 10 30\"\nDEFAULT: 5 50\n", 
                    "id": "#macs2-callpeak.cwl/mfold"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "For staging in a case of trigger is set to false", 
                    "id": "#macs2-callpeak.cwl/moder_r_file_staged"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 999, 
                        "prefix": "-n", 
                        "valueFrom": "${\n    if (self == null || inputs.trigger == false){\n      return default_name();\n    } else {\n      return self;\n    }\n}\n"
                    }, 
                    "default": null, 
                    "doc": "The name string of the experiment. MACS will use this string NAME to create output files like \u2018NAME_peaks.xls\u2019,\n\u2018NAME_negative_peaks.xls\u2019, \u2018NAME_peaks.bed\u2019 , \u2018NAME_summits.bed\u2019, \u2018NAME_model.r\u2019 and so on.\nSo please avoid any confliction between these filenames and your existing files.\nDEFAULT: generated on the base of the treatment input\n", 
                    "id": "#macs2-callpeak.cwl/name"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "For staging in a case of trigger is set to false", 
                    "id": "#macs2-callpeak.cwl/narrow_peak_file_staged"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 33, 
                        "prefix": "--nolambda"
                    }, 
                    "doc": "If True, MACS will use fixed background lambda as local lambda for every\npeak region. Normally, MACS calculates a dynamic local lambda to reflect the\nlocal bias due to potential chromatin structure.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/nolambda"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 24, 
                        "prefix": "--nomodel"
                    }, 
                    "doc": "Whether or not to build the shifting model. If True,  MACS will not build\nmodel. by default it means  shifting size = 100, try to set extsize to change it.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/nomodel"
                }, 
                {
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "inputBinding": {
                        "position": 28, 
                        "prefix": "-p"
                    }, 
                    "doc": "Pvalue cutoff for peak detection. DEFAULT: not set.  -q, and -p are mutually\nexclusive. If pvalue cutoff is set, qvalue will not be calculated and reported\nas -1  in the final .xls file.\nDEFAULT: null\n", 
                    "id": "#macs2-callpeak.cwl/p_value"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "For staging in a case of trigger is set to false", 
                    "id": "#macs2-callpeak.cwl/peak_summits_file_staged"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "For staging in a case of trigger is set to false", 
                    "id": "#macs2-callpeak.cwl/peak_xls_file_staged"
                }, 
                {
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "inputBinding": {
                        "position": 27, 
                        "prefix": "-q"
                    }, 
                    "doc": "Minimum FDR (q-value) cutoff for peak detection. -q, and\n-p are mutually exclusive.\nDEFAULT: 0.05\n", 
                    "id": "#macs2-callpeak.cwl/q_value"
                }, 
                {
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "inputBinding": {
                        "position": 30, 
                        "prefix": "--ratio"
                    }, 
                    "doc": "When set, use a custom scaling ratio of ChIP/control\n(e.g. calculated using NCIS) for linear scaling.\nDEFAULT: null\n", 
                    "id": "#macs2-callpeak.cwl/ratio"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "default": "#!/bin/bash\nif [ \"$0\" = True ]\nthen\n  echo \"Run: macs2 callpeak \" ${@:1}\n  ls | grep -v ${@: -1}.log | xargs rm\n  macs2 callpeak \"${@:1}\"\nelse\n  echo \"Skip macs2 callpeak \" ${@:1}\nfi\n", 
                    "inputBinding": {
                        "position": 1
                    }, 
                    "doc": "Bash function to run MACS2 callpeak with all input parameters or skip it if trigger is false\n", 
                    "id": "#macs2-callpeak.cwl/script"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 32, 
                        "prefix": "--seed"
                    }, 
                    "doc": "Set the random seed while down sampling data. Must be\na non-negative integer in order to be effective.\nDEFAULT: null\n", 
                    "id": "#macs2-callpeak.cwl/seed"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 25, 
                        "prefix": "--shift"
                    }, 
                    "doc": "(NOT the legacy --shiftsize option!) The arbitrary shift in bp. Use discretion\nwhile setting it other than default value. When NOMODEL is set, MACS will use\nthis value to move cutting ends (5'') towards 5''->3'' direction then apply\nEXTSIZE to extend them to fragments. When this value is negative, ends will\nbe moved toward 3''->5'' direction. Recommended to keep it as default 0 for\nChIP-Seq datasets, or -1 * half of EXTSIZE together with EXTSIZE option for\ndetecting enriched cutting loci such as certain DNAseI-Seq datasets. Note, you\ncan''t set values other than 0 if format is BAMPE for paired-end data.\nDEFAULT: 0\n", 
                    "id": "#macs2-callpeak.cwl/shift"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 34, 
                        "prefix": "--slocal"
                    }, 
                    "doc": "The small nearby region in basepairs to calculate\ndynamic lambda. This is used to capture the bias near\nthe peak summit region. Invalid if there is no control\ndata. If you set this to 0, MACS will skip slocal\nlambda calculation. *Note* that MACS will always\nperform a d-size local lambda calculation. The final\nlocal bias should be the maximum of the lambda value\nfrom d, slocal, and llocal size windows.\nDEFAULT: 1000\n", 
                    "id": "#macs2-callpeak.cwl/slocal"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 19, 
                        "prefix": "--SPMR"
                    }, 
                    "doc": "If True, MACS will save signal per million reads for fragment pileup profiles.\nRequire --bdg to be set.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/spmr"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 29, 
                        "prefix": "--to-large"
                    }, 
                    "doc": "When set, scale the small sample up to the bigger sample. By default, the\nbigger dataset will be scaled down towards the smaller dataset, which will lead\nto smaller p/qvalues and more specific results. Keep in mind that scaling down\nwill bring down background noise more.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/to_large"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 18, 
                        "prefix": "--trackline"
                    }, 
                    "doc": "Tells MACS to include trackline with bedGraph files. To include this trackline\nwhile displaying bedGraph at UCSC genome browser, can show name and description\nof the file as well. Require -B to be set.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/trackline"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "doc": "For staging in a case of trigger is set to false", 
                    "id": "#macs2-callpeak.cwl/treat_pileup_bdg_file_staged"
                }, 
                {
                    "type": [
                        "File", 
                        {
                            "type": "array", 
                            "items": "File"
                        }
                    ], 
                    "inputBinding": {
                        "position": 10, 
                        "prefix": "-t"
                    }, 
                    "doc": "This is the only REQUIRED parameter for MACS. File can be in any supported format specified by \u2013format option.\nCheck \u2013format for detail. If you have more than one alignment files, you can specify them as `-t A B C`.\nMACS will pool up all these files together.\n", 
                    "id": "#macs2-callpeak.cwl/treatment"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "default": true, 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "doc": "If true - run MACS2, if false - return staged files\n", 
                    "id": "#macs2-callpeak.cwl/trigger"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 20, 
                        "prefix": "--tsize"
                    }, 
                    "doc": "Tag size. This will overide the auto detected tag size.\nDEFAULT: False\n", 
                    "id": "#macs2-callpeak.cwl/tsize"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 41, 
                        "prefix": "--verbose"
                    }, 
                    "doc": "Log level\n", 
                    "id": "#macs2-callpeak.cwl/verbose"
                }
            ], 
            "outputs": [
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${ if (inputs.name == null || inputs.trigger == false){ return default_name(inputs.broad_peak_file_staged, '_peaks.broadPeak'); } else { return inputs.name + '_peaks.broadPeak'; } }"
                    }, 
                    "id": "#macs2-callpeak.cwl/broad_peak_file"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${ if (inputs.name == null || inputs.trigger == false){ return default_name(inputs.control_lambda_bdg_file_staged, '_control_lambda.bdg'); } else { return inputs.name + '_control_lambda.bdg'; } }"
                    }, 
                    "id": "#macs2-callpeak.cwl/control_lambda_bdg_file"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${ if (inputs.name == null || inputs.trigger == false){ return default_name(inputs.gapped_peak_file_staged, '_peaks.gappedPeak'); } else { return inputs.name + '_peaks.gappedPeak'; } }"
                    }, 
                    "id": "#macs2-callpeak.cwl/gapped_peak_file"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${\n    if (inputs.name == null || inputs.trigger == false){\n      return default_name(null, '.log');\n    } else {\n      return inputs.name + '.log';\n    }\n}\n"
                    }, 
                    "id": "#macs2-callpeak.cwl/macs_log"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${ if (inputs.name == null || inputs.trigger == false){ return default_name(inputs.moder_r_file_staged, '_model.r'); } else { return inputs.name + '_model.r'; } }"
                    }, 
                    "id": "#macs2-callpeak.cwl/moder_r_file"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${ if (inputs.name == null || inputs.trigger == false){ return default_name(inputs.narrow_peak_file_staged, '_peaks.narrowPeak'); } else { return inputs.name + '_peaks.narrowPeak'; } }"
                    }, 
                    "id": "#macs2-callpeak.cwl/narrow_peak_file"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${ if (inputs.name == null || inputs.trigger == false){ return default_name(inputs.peak_summits_file_staged, '_summits.bed'); } else { return inputs.name + '_summits.bed'; } }"
                    }, 
                    "id": "#macs2-callpeak.cwl/peak_summits_file"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${ if (inputs.name == null || inputs.trigger == false){ return default_name(inputs.peak_xls_file_staged, '_peaks.xls'); } else { return inputs.name + '_peaks.xls'; } }"
                    }, 
                    "id": "#macs2-callpeak.cwl/peak_xls_file"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "outputBinding": {
                        "glob": "${ if (inputs.name == null || inputs.trigger == false){ return default_name(inputs.treat_pileup_bdg_file_staged, '_treat_pileup.bdg'); } else { return inputs.name + '_treat_pileup.bdg'; } }"
                    }, 
                    "id": "#macs2-callpeak.cwl/treat_pileup_bdg_file"
                }
            ], 
            "baseCommand": [
                "bash", 
                "-c"
            ], 
            "arguments": [
                {
                    "valueFrom": "${ if (inputs.name == null || inputs.trigger == false ){ return ' 2> ' + default_name(null, '.log'); } else { return ' 2> ' + inputs.name + '.log'; } }", 
                    "position": 100000, 
                    "shellQuote": false
                }
            ], 
            "s:name": "macs2-callpeak", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "Tool is used to perform peak calling using MACS2.\nInput Trigger (default: true) allows to skip all calculation and return\nall input files unchanged. To set files to be returned in case of Trigger == false,\nuse the following inputs:\n  peak_xls_file_staged:\n  narrow_peak_file_staged:\n  broad_peak_file_staged:\n  gapped_peak_file_staged:\n  peak_summits_file_staged:\n  moder_r_file_staged:\n  treat_pileup_bdg_file_staged:\n  control_lambda_bdg_file_staged:\n", 
            "id": "#macs2-callpeak.cwl", 
            "http://schema.org/mainEntity": {
                "class": "http://schema.org/SoftwareSourceCode", 
                "s:about": "With the improvement of sequencing techniques, chromatin immunoprecipitation followed by high throughput sequencing (ChIP-Seq) is getting popular to study genome-wide protein-DNA interactions. To address the lack of powerful ChIP-Seq analysis method, we present a novel algorithm, named Model-based Analysis of ChIP-Seq (MACS), for identifying transcript factor binding sites. MACS captures the influence of genome complexity to evaluate the significance of enriched ChIP regions, and MACS improves the spatial resolution of binding sites through combining the information of both sequencing tag position and orientation. MACS can be easily used for ChIP-Seq data alone, or with control sample with the increase of specificity.\n", 
                "s:codeRepository": "https://github.com/taoliu/MACS", 
                "s:targetProduct": {
                    "class": "http://schema.org/SoftwareApplication", 
                    "s:applicationCategory": "command line tool", 
                    "http://schema.org/softwareVersion": "2.1.1.20160309"
                }, 
                "s:publication": [
                    {
                        "class": "http://schema.org/ScholarlyArticle", 
                        "id": "#macs2-metadata.yaml/10.1186/gb-2008-9-9-r137", 
                        "s:author": [
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Zhang Y"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Liu T"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Meyer CA"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Eeckhoute J"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Johnson DS"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Bernstein BE"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Nusbaum C"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Myers RM"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Brown M"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Li W"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "http://schema.org/name": "Liu XS"
                            }
                        ], 
                        "datePublished": "17 September 2008", 
                        "http://schema.org/name": "Model-based analysis of ChIP-Seq (MACS)", 
                        "http://schema.org/url": "https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137"
                    }
                ], 
                "s:creator": [
                    {
                        "class": "http://schema.org/Organization", 
                        "s:member": [
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "Main author", 
                                "http://schema.org/name": "Yong Zhang"
                            }, 
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "Main author", 
                                "http://schema.org/name": "Tao Liu"
                            }
                        ], 
                        "http://schema.org/name": "Xiaole Shirley Liu's Lab"
                    }
                ], 
                "id": "#macs2-metadata.yaml", 
                "name": "#macs2-metadata.yaml", 
                "http://schema.org/name": "MACS2", 
                "http://schema.org/url": "http://liulab.dfci.harvard.edu/MACS/", 
                "http://schema.org/license": [
                    "https://opensource.org/licenses/BSD-3-Clause"
                ], 
                "http://schema.org/programmingLanguage": "Python", 
                "http://schema.org/discussionUrl": [
                    "https://groups.google.com/forum/#!forum/macs-announcement", 
                    "https://github.com/taoliu/MACS/issues"
                ]
            }, 
            "http://schema.org/downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/macs2-callpeak.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:michael.kotliar@cchmc.org", 
                                            "http://schema.org/name": "Michael Kotliar", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0002-6486-3898"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "usage: macs2 callpeak [-h] -t TFILE [TFILE ...] [-c [CFILE [CFILE ...]]]\n                      [-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}]\n                      [-g GSIZE] [--keep-dup KEEPDUPLICATES]\n                      [--buffer-size BUFFER_SIZE] [--outdir OUTDIR] [-n NAME]\n                      [-B] [--verbose VERBOSE] [--trackline] [--SPMR]\n                      [-s TSIZE] [--bw BW] [-m MFOLD MFOLD] [--fix-bimodal]\n                      [--nomodel] [--shift SHIFT] [--extsize EXTSIZE]\n                      [-q QVALUE | -p PVALUE] [--to-large] [--ratio RATIO]\n                      [--down-sample] [--seed SEED] [--tempdir TEMPDIR]\n                      [--nolambda] [--slocal SMALLLOCAL] [--llocal LARGELOCAL]\n                      [--broad] [--broad-cutoff BROADCUTOFF]\n                      [--cutoff-analysis] [--call-summits]\n                      [--fe-cutoff FECUTOFF]\n\noptional arguments:\n  -h, --help            show this help message and exit\n\nInput files arguments:\n  -t TFILE [TFILE ...], --treatment TFILE [TFILE ...]\n                        ChIP-seq treatment file. If multiple files are given\n                        as '-t A B C', then they will all be read and pooled\n                        together. REQUIRED.\n  -c [CFILE [CFILE ...]], --control [CFILE [CFILE ...]]\n                        Control file. If multiple files are given as '-c A B\n                        C', they will be pooled to estimate ChIP-seq\n                        background noise.\n  -f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}, --format {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}\n                        Format of tag file, \"AUTO\", \"BED\" or \"ELAND\" or\n                        \"ELANDMULTI\" or \"ELANDEXPORT\" or \"SAM\" or \"BAM\" or\n                        \"BOWTIE\" or \"BAMPE\" or \"BEDPE\". The default AUTO\n                        option will let MACS decide which format (except for\n                        BAMPE and BEDPE which should be implicitly set) the\n                        file is. Please check the definition in README. Please\n                        note that if the format is set as BAMPE or BEDPE,\n                        MACS2 will call its special Paired-end mode to call\n                        peaks by piling up the actual ChIPed fragments defined\n                        by both aligned ends, instead of predicting the\n                        fragment size first and extending reads. Also please\n                        note that the BEDPE only contains three columns, and\n                        is NOT the same BEDPE format used by BEDTOOLS.\n                        DEFAULT: \"AUTO\"\n  -g GSIZE, --gsize GSIZE\n                        Effective genome size. It can be 1.0e+9 or 1000000000,\n                        or shortcuts:'hs' for human (2.7e9), 'mm' for mouse\n                        (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for\n                        fruitfly (1.2e8), Default:hs\n  --keep-dup KEEPDUPLICATES\n                        It controls the MACS behavior towards duplicate tags\n                        at the exact same location -- the same coordination\n                        and the same strand. The 'auto' option makes MACS\n                        calculate the maximum tags at the exact same location\n                        based on binomal distribution using 1e-5 as pvalue\n                        cutoff; and the 'all' option keeps every tags. If an\n                        integer is given, at most this number of tags will be\n                        kept at the same location. Note, if you've used\n                        samtools or picard to flag reads as 'PCR/Optical\n                        duplicate' in bit 1024, MACS2 will still read them\n                        although the reads may be decided by MACS2 as\n                        duplicate later. The default is to keep one tag at the\n                        same location. Default: 1\n  --buffer-size BUFFER_SIZE\n                        Buffer size for incrementally increasing internal\n                        array size to store reads alignment information. In\n                        most cases, you don't have to change this parameter.\n                        However, if there are large number of\n                        chromosomes/contigs/scaffolds in your alignment, it's\n                        recommended to specify a smaller buffer size in order\n                        to decrease memory usage (but it will take longer time\n                        to read alignment files). Minimum memory requested for\n                        reading an alignment file is about # of CHROMOSOME *\n                        BUFFER_SIZE * 2 Bytes. DEFAULT: 100000\n\nOutput arguments:\n  --outdir OUTDIR       If specified all output files will be written to that\n                        directory. Default: the current working directory\n  -n NAME, --name NAME  Experiment name, which will be used to generate output\n                        file names. DEFAULT: \"NA\"\n  -B, --bdg             Whether or not to save extended fragment pileup, and\n                        local lambda tracks (two files) at every bp into a\n                        bedGraph file. DEFAULT: False\n  --verbose VERBOSE     Set verbose level of runtime message. 0: only show\n                        critical message, 1: show additional warning message,\n                        2: show process information, 3: show debug messages.\n                        DEFAULT:2\n  --trackline           Tells MACS to include trackline with bedGraph files.\n                        To include this trackline while displaying bedGraph at\n                        UCSC genome browser, can show name and description of\n                        the file as well. However my suggestion is to convert\n                        bedGraph to bigWig, then show the smaller and faster\n                        binary bigWig file at UCSC genome browser, as well as\n                        downstream analysis. Require -B to be set. Default:\n                        Not include trackline.\n  --SPMR                If True, MACS will save signal per million reads for\n                        fragment pileup profiles. Require -B to be set.\n                        Default: False\n\nShifting model arguments:\n  -s TSIZE, --tsize TSIZE\n                        Tag size. This will override the auto detected tag\n                        size. DEFAULT: Not set\n  --bw BW               Band width for picking regions to compute fragment\n                        size. This value is only used while building the\n                        shifting model. DEFAULT: 300\n  -m MFOLD MFOLD, --mfold MFOLD MFOLD\n                        Select the regions within MFOLD range of high-\n                        confidence enrichment ratio against background to\n                        build model. Fold-enrichment in regions must be lower\n                        than upper limit, and higher than the lower limit. Use\n                        as \"-m 10 30\". DEFAULT:5 50\n  --fix-bimodal         Whether turn on the auto pair model process. If set,\n                        when MACS failed to build paired model, it will use\n                        the nomodel settings, the --exsize parameter to extend\n                        each tags towards 3' direction. Not to use this\n                        automate fixation is a default behavior now. DEFAULT:\n                        False\n  --nomodel             Whether or not to build the shifting model. If True,\n                        MACS will not build model. by default it means\n                        shifting size = 100, try to set extsize to change it.\n                        DEFAULT: False\n  --shift SHIFT         (NOT the legacy --shiftsize option!) The arbitrary\n                        shift in bp. Use discretion while setting it other\n                        than default value. When NOMODEL is set, MACS will use\n                        this value to move cutting ends (5') towards 5'->3'\n                        direction then apply EXTSIZE to extend them to\n                        fragments. When this value is negative, ends will be\n                        moved toward 3'->5' direction. Recommended to keep it\n                        as default 0 for ChIP-Seq datasets, or -1 * half of\n                        EXTSIZE together with EXTSIZE option for detecting\n                        enriched cutting loci such as certain DNAseI-Seq\n                        datasets. Note, you can't set values other than 0 if\n                        format is BAMPE or BEDPE for paired-end data. DEFAULT:\n                        0.\n  --extsize EXTSIZE     The arbitrary extension size in bp. When nomodel is\n                        true, MACS will use this value as fragment size to\n                        extend each read towards 3' end, then pile them up.\n                        It's exactly twice the number of obsolete SHIFTSIZE.\n                        In previous language, each read is moved 5'->3'\n                        direction to middle of fragment by 1/2 d, then\n                        extended to both direction with 1/2 d. This is\n                        equivalent to say each read is extended towards 5'->3'\n                        into a d size fragment. DEFAULT: 200. EXTSIZE and\n                        SHIFT can be combined when necessary. Check SHIFT\n                        option.\n\nPeak calling arguments:\n  -q QVALUE, --qvalue QVALUE\n                        Minimum FDR (q-value) cutoff for peak detection.\n                        DEFAULT: 0.05. -q, and -p are mutually exclusive.\n  -p PVALUE, --pvalue PVALUE\n                        Pvalue cutoff for peak detection. DEFAULT: not set.\n                        -q, and -p are mutually exclusive. If pvalue cutoff is\n                        set, qvalue will not be calculated and reported as -1\n                        in the final .xls file.\n  --to-large            When set, scale the small sample up to the bigger\n                        sample. By default, the bigger dataset will be scaled\n                        down towards the smaller dataset, which will lead to\n                        smaller p/qvalues and more specific results. Keep in\n                        mind that scaling down will bring down background\n                        noise more. DEFAULT: False\n  --ratio RATIO         When set, use a custom scaling ratio of ChIP/control\n                        (e.g. calculated using NCIS) for linear scaling.\n                        DEFAULT: ingore\n  --down-sample         When set, random sampling method will scale down the\n                        bigger sample. By default, MACS uses linear scaling.\n                        Warning: This option will make your result unstable\n                        and irreproducible since each time, random reads would\n                        be selected. Consider to use 'randsample' script\n                        instead. <not implmented>If used together with --SPMR,\n                        1 million unique reads will be randomly picked.</not\n                        implemented> Caution: due to the implementation, the\n                        final number of selected reads may not be as you\n                        expected! DEFAULT: False\n  --seed SEED           Set the random seed while down sampling data. Must be\n                        a non-negative integer in order to be effective.\n                        DEFAULT: not set\n  --tempdir TEMPDIR     Optional directory to store temp files. DEFAULT: /tmp\n  --nolambda            If True, MACS will use fixed background lambda as\n                        local lambda for every peak region. Normally, MACS\n                        calculates a dynamic local lambda to reflect the local\n                        bias due to potential chromatin structure.\n  --slocal SMALLLOCAL   The small nearby region in basepairs to calculate\n                        dynamic lambda. This is used to capture the bias near\n                        the peak summit region. Invalid if there is no control\n                        data. If you set this to 0, MACS will skip slocal\n                        lambda calculation. *Note* that MACS will always\n                        perform a d-size local lambda calculation. The final\n                        local bias should be the maximum of the lambda value\n                        from d, slocal, and llocal size windows. DEFAULT: 1000\n  --llocal LARGELOCAL   The large nearby region in basepairs to calculate\n                        dynamic lambda. This is used to capture the surround\n                        bias. If you set this to 0, MACS will skip llocal\n                        lambda calculation. *Note* that MACS will always\n                        perform a d-size local lambda calculation. The final\n                        local bias should be the maximum of the lambda value\n                        from d, slocal, and llocal size windows. DEFAULT:\n                        10000.\n  --broad               If set, MACS will try to call broad peaks by linking\n                        nearby highly enriched regions. The linking region is\n                        controlled by another cutoff through --linking-cutoff.\n                        The maximum linking region length is 4 times of d from\n                        MACS. DEFAULT: False\n  --broad-cutoff BROADCUTOFF\n                        Cutoff for broad region. This option is not available\n                        unless --broad is set. If -p is set, this is a pvalue\n                        cutoff, otherwise, it's a qvalue cutoff. DEFAULT: 0.1\n  --cutoff-analysis     While set, MACS2 will analyze number or total length\n                        of peaks that can be called by different p-value\n                        cutoff then output a summary table to help user decide\n                        a better cutoff. The table will be saved in\n                        NAME_cutoff_analysis.txt file. Note, minlen and maxgap\n                        may affect the results. WARNING: May take ~30 folds\n                        longer time to finish. DEFAULT: False\n\nPost-processing options:\n  --call-summits        If set, MACS will use a more sophisticated signal\n                        processing approach to find subpeak summits in each\n                        enriched peak region. DEFAULT: False\n  --fe-cutoff FECUTOFF  When set, the value will be used to filter out peaks\n                        with low fold-enrichment. Note, MACS2 use 1.0 as\n                        pseudocount while calculating fold-enrichment.\n                        DEFAULT: 1.0\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "ShellCommandRequirement"
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/scidap:v0.0.2", 
                    "dockerFile": "$import: ../dockerfiles/scidap/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 6
                    }, 
                    "doc": "Output file from MACS2 peak calling\n", 
                    "id": "#macs2-island-count.cwl/input_file"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "default": "#!/usr/bin/env python\nimport sys, re\nfragments, islands = 0, 0\nwith open(sys.argv[1], 'r') as infile:\n  for line in infile:\n      if re.match('^# d = ', line):\n          fragments = int(line.split('d = ')[1])\n          continue\n      if re.match('^#', line):\n          continue\n      if line.strip() != \"\":\n          islands = islands + 1\nislands = islands - 1\nprint fragments, '\\n', islands\n", 
                    "inputBinding": {
                        "position": 5
                    }, 
                    "doc": "Python script to get ISLANDS and FRAGMENTS from MACS2 output\n", 
                    "id": "#macs2-island-count.cwl/script"
                }
            ], 
            "outputs": [
                {
                    "type": "int", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "island_count.log", 
                        "outputEval": "$(parseInt(self[0].contents.split(/\\r?\\n/)[0]))"
                    }, 
                    "id": "#macs2-island-count.cwl/fragments"
                }, 
                {
                    "type": "int", 
                    "outputBinding": {
                        "loadContents": true, 
                        "glob": "island_count.log", 
                        "outputEval": "$(parseInt(self[0].contents.split(/\\r?\\n/)[1]))"
                    }, 
                    "id": "#macs2-island-count.cwl/islands"
                }
            ], 
            "baseCommand": [
                "python", 
                "-c"
            ], 
            "arguments": [
                {
                    "valueFrom": "$(\" > island_count.log\")", 
                    "position": 100000, 
                    "shellQuote": false
                }
            ], 
            "s:name": "macs2-island-count", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "Tool is used to return an estimated fragment size and islands count from\nxls file generated by MACS2 callpeak\n", 
            "id": "#macs2-island-count.cwl", 
            "http://schema.org/mainEntity": {
                "$import": "#macs2-metadata.yaml"
            }, 
            "http://schema.org/downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/macs2-island-count.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:michael.kotliar@cchmc.org", 
                                            "http://schema.org/name": "Michael Kotliar", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0002-6486-3898"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "Runs python code from the script input\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "ShellCommandRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement", 
                    "expressionLib": [
                        "var default_output_filename = function() { return inputs.bowtie_log.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".stat\"; };"
                    ]
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/scidap:v0.0.2", 
                    "dockerFile": "$import: ../dockerfiles/scidap/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 6
                    }, 
                    "doc": "Log file from Bowtie\n", 
                    "id": "#python-get-stat.cwl/bowtie_log"
                }, 
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 7
                    }, 
                    "doc": "Log file from samtools rmdup\n", 
                    "id": "#python-get-stat.cwl/rmdup_log"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "default": "#!/usr/bin/env python\nimport sys, re\nTOTAL, ALIGNED, SUPRESSED, USED = 100, 80, 0, 0\nwith open(sys.argv[1], 'r') as bowtie_log:\n  for line in bowtie_log:\n    if 'processed:' in line:\n      TOTAL = int(line.split('processed:')[1])\n    if 'alignment:' in line:\n      ALIGNED = int(line.split('alignment:')[1].split()[0])\n    if 'due to -m:' in line:\n      SUPRESSED = int(line.split('due to -m:')[1].split()[0])\nUSED = ALIGNED\nwith open(sys.argv[2], 'r') as rmdup_log:\n  for line in rmdup_log:\n    if '/' in line and 'Skip' not in line:\n      splt = line.split('/')\n      USED = int((splt[1].split('='))[0].strip()) - int((splt[0].split(']'))[1].strip())\nprint TOTAL, ALIGNED, SUPRESSED, USED\n", 
                    "inputBinding": {
                        "position": 5
                    }, 
                    "doc": "Python script to get TOTAL, ALIGNED, SUPRESSED, USED values from log files\n", 
                    "id": "#python-get-stat.cwl/script"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "$(default_output_filename())"
                    }, 
                    "id": "#python-get-stat.cwl/output"
                }
            ], 
            "baseCommand": [
                "python", 
                "-c"
            ], 
            "arguments": [
                {
                    "valueFrom": "$(\" > \" + default_output_filename())", 
                    "position": 100000, 
                    "shellQuote": false
                }
            ], 
            "s:downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/python-get-stat.cwl", 
            "s:license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "s:creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:michael.kotliar@cchmc.org", 
                                            "http://schema.org/name": "Michael Kotliar", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0002-6486-3898"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "doc": "Tool to process and combine log files generated by Bowtie aligner and samtools rmdup.\n", 
            "id": "#python-get-stat.cwl", 
            "http://schema.org/name": "python-get-stat", 
            "http://schema.org/codeRepository": "https://github.com/SciDAP/workflows", 
            "http://schema.org/isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "http://schema.org/about": "Runs python code from the script input\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "ShellCommandRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement", 
                    "expressionLib": [
                        "var default_output_filename = function() { if (inputs.trigger == true){ return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".r.bam\"; } else { return inputs.input_file.location.split('/').slice(-1)[0]; } };"
                    ]
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.input_file)"
                    ]
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/samtools:v1.4", 
                    "dockerFile": "$import: ../dockerfiles/samtools/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "default": "#!/bin/bash\nif [ \"$0\" = True ]\nthen\n  echo \"Run: samtools rmdup \" ${@:1}\n  samtools rmdup \"${@:1}\"\nelse\n  echo \"Skip samtools rmdup \" ${@:1}\nfi\n", 
                    "inputBinding": {
                        "position": 1
                    }, 
                    "doc": "Bash function to run samtools rmdup with all input parameters or skip it if trigger is false\n", 
                    "id": "#samtools-rmdup.cwl/bash_script"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 4, 
                        "prefix": "-S"
                    }, 
                    "doc": "treat PE reads as SE in rmdup (force -s)\n", 
                    "id": "#samtools-rmdup.cwl/force_single_end"
                }, 
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 10
                    }, 
                    "doc": "Input sorted bam file.\n", 
                    "id": "#samtools-rmdup.cwl/input_file"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 11, 
                        "valueFrom": "${\n    if (self == null || inputs.trigger == false){\n      return default_output_filename();\n    } else {\n      return self;\n    }\n}\n"
                    }, 
                    "default": null, 
                    "doc": "Writes the output bam file to output_filename if set,\notherwise generates output_filename on the base of input_file\n", 
                    "id": "#samtools-rmdup.cwl/output_filename"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 3, 
                        "prefix": "-s"
                    }, 
                    "doc": "rmdup for SE reads\n", 
                    "id": "#samtools-rmdup.cwl/single_end"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "default": true, 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "doc": "If true - run samtools rmdup, if false - return input_file, previously staged into output directory\n", 
                    "id": "#samtools-rmdup.cwl/trigger"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "${\n  if (inputs.output_filename == null || inputs.trigger == false){\n    return default_output_filename() + '.rmdup';\n  } else {\n    return inputs.output_filename + '.rmdup';\n  }\n}\n"
                    }, 
                    "id": "#samtools-rmdup.cwl/rmdup_log"
                }, 
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "${\n    if (inputs.output_filename == null || inputs.trigger == false){\n      return default_output_filename();\n    } else {\n      return inputs.output_filename;\n    }\n}\n"
                    }, 
                    "secondaryFiles": "${\n    if (inputs.input_file.secondaryFiles && inputs.trigger == false){\n      return inputs.input_file.secondaryFiles;\n    } else {\n      return \"null\";\n    }\n  }\n", 
                    "doc": "File with removed duplicates or input_file with optional secondaryFiles", 
                    "id": "#samtools-rmdup.cwl/rmdup_output"
                }
            ], 
            "baseCommand": [
                "bash", 
                "-c"
            ], 
            "arguments": [
                {
                    "valueFrom": "${\n  if (inputs.output_filename == null || inputs.trigger == false){\n    return \" > \" + default_output_filename() + \".rmdup 2>&1\";\n  } else {\n    return \" > \" + inputs.output_filename + \".rmdup 2>&1\";\n  }\n}\n", 
                    "position": 100000, 
                    "shellQuote": false
                }
            ], 
            "s:name": "samtools-rmdup", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "This tool is used to remove duplicates from input coordinate sorted BAM file\nInput Trigger (default: true) allows to skip all calculation and return\nall input files unchanged. To set files to be returned in case of Trigger == false,\nuse the following inputs:\n  input_file\n", 
            "id": "#samtools-rmdup.cwl", 
            "http://schema.org/mainEntity": {
                "class": "http://schema.org/SoftwareSourceCode", 
                "s:about": "A suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories: Samtools (Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format), BCFtools (Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants) and HTSlib (A C library for reading/writing high-throughput sequencing data).\n", 
                "s:codeRepository": "https://github.com/samtools/samtools.git", 
                "s:targetProduct": {
                    "class": "http://schema.org/SoftwareApplication", 
                    "s:applicationCategory": "commandline tool", 
                    "http://schema.org/softwareVersion": "1.4"
                }, 
                "s:publication": [
                    {
                        "class": "http://schema.org/ScholarlyArticle", 
                        "id": "#btr509", 
                        "s:url": "http://www.ncbi.nlm.nih.gov/pubmed/21903627", 
                        "http://schema.org/name": "(Li, 2011) A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics."
                    }, 
                    {
                        "class": "http://schema.org/ScholarlyArticle", 
                        "id": "#btr076", 
                        "s:url": "http://www.ncbi.nlm.nih.gov/pubmed/21320865", 
                        "http://schema.org/name": "(Li, 2011) Improving SNP discovery by base alignment quality. Bioinformatics."
                    }, 
                    {
                        "class": "http://schema.org/ScholarlyArticle", 
                        "id": "#btp352", 
                        "s:url": "http://www.ncbi.nlm.nih.gov/pubmed/19505943", 
                        "http://schema.org/name": "(Li et al., 2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics."
                    }
                ], 
                "s:creator": [
                    {
                        "class": "http://schema.org/Organization", 
                        "s:member": [
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "wrote most of the initial source codes of SAMtools and various converters.", 
                                "http://schema.org/name": "Heng Li"
                            }
                        ], 
                        "http://schema.org/name": "Sanger Institute"
                    }, 
                    {
                        "class": "http://schema.org/Organization", 
                        "s:member": [
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "A major contributor to the\nSAM/BAM specification. He designed and implemented the BGZF format, the\nunderlying indexable compression format for the BAM format. BGZF does\nnot support arithmetic between file offsets.\n", 
                                "http://schema.org/name": "Bob Handsaker"
                            }
                        ], 
                        "http://schema.org/name": "Broad Institute"
                    }, 
                    {
                        "class": "http://schema.org/Organization", 
                        "s:member": [
                            {
                                "class": "http://schema.org/Person", 
                                "s:description": "Designed and implemented the\nRAZF format, an alternative indexable compression format. RAZF is no longer\nused by or provided with SAMtools. Source code remains available in older\nSAMtools 0.1.x releases and from the standalone branch in the repository.\n", 
                                "http://schema.org/name": "Jue Ruan"
                            }
                        ], 
                        "http://schema.org/name": "Beijing Genome Institute"
                    }, 
                    {
                        "class": "http://schema.org/Person", 
                        "s:description": "updated novo2sam.pl to support gapped alignment by novoalign.", 
                        "http://schema.org/name": "Colin Hercus"
                    }, 
                    {
                        "class": "http://schema.org/Person", 
                        "s:description": "contributed the header parsing library sam_header.c and sam2vcf.pl script.", 
                        "http://schema.org/name": "Petr Danecek"
                    }
                ], 
                "id": "#samtools-metadata.yaml", 
                "name": "#samtools-metadata.yaml", 
                "http://schema.org/name": "samtools", 
                "http://schema.org/url": "http://www.htslib.org/", 
                "http://schema.org/license": [
                    "https://opensource.org/licenses/MIT", 
                    "https://opensource.org/licenses/BSD-3-Clause"
                ], 
                "http://schema.org/programmingLanguage": "C, Perl", 
                "http://schema.org/discussionUrl": [
                    "https://lists.sourceforge.net/lists/listinfo/samtools-help", 
                    "https://lists.sourceforge.net/lists/listinfo/samtools-devel"
                ]
            }, 
            "http://schema.org/downloadUrl": "https://github.com/SciDAP/workflows/blob/master/tools/samtools-rmdup.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:michael.kotliar@cchmc.org", 
                                            "http://schema.org/name": "Michael Kotliar", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0002-6486-3898"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "Usage:  samtools rmdup [-sS] <input.srt.bam> <output.bam>\nOption: -s    rmdup for SE reads\n        -S    treat PE reads as SE in rmdup (force -s)\n      --input-fmt-option OPT[=VAL]\n               Specify a single input file format option in the form\n               of OPTION or OPTION=VALUE\n      --output-fmt FORMAT[,OPT[=VAL]]...\n               Specify output format (SAM, BAM, CRAM)\n      --output-fmt-option OPT[=VAL]\n               Specify a single output file format option in the form\n               of OPTION or OPTION=VALUE\n      --reference FILE\n               Reference sequence FASTA FILE [null]\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "$import": "#envvar-global.yml"
                }, 
                {
                    "class": "ShellCommandRequirement"
                }, 
                {
                    "class": "InitialWorkDirRequirement", 
                    "listing": [
                        "$(inputs.sort_input)"
                    ]
                }, 
                {
                    "class": "InlineJavascriptRequirement", 
                    "expressionLib": [
                        "var ext = function() { if (inputs.csi && !inputs.bai){ return '.csi'; } else { return '.bai'; } };", 
                        "var default_bam = function() { if (inputs.trigger == true){ return inputs.sort_input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".sorted.bam\"; } else { return inputs.sort_input.location.split('/').slice(-1)[0]; } };"
                    ]
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "biowardrobe2/samtools:v1.4", 
                    "dockerFile": "$import: ../dockerfiles/samtools/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "Generate BAI-format index for BAM files [default]. If input isn't cram.\n", 
                    "id": "#samtools-sort-index.cwl/bai"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "default": "#!/bin/bash\nif [ \"$0\" = True ]\nthen\n  echo \"Run: samtools index \" ${@:1}\n  samtools index \"${@:1}\"\nelse\n  echo \"Skip samtools index \" ${@:1}\nfi\n", 
                    "inputBinding": {
                        "position": 20
                    }, 
                    "doc": "Bash function to run samtools index with all input parameters or skip it if trigger is false\n", 
                    "id": "#samtools-sort-index.cwl/bash_script_index"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "default": "#!/bin/bash\nif [ \"$0\" = True ]\nthen\n  echo \"Run: samtools sort \" ${@:1}\n  samtools sort \"${@:1}\"\nelse\n  echo \"Skip samtools sort \" ${@:1}\nfi\n", 
                    "inputBinding": {
                        "position": 5
                    }, 
                    "doc": "Bash function to run samtools sort with all input parameters or skip it if trigger is false\n", 
                    "id": "#samtools-sort-index.cwl/bash_script_sort"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "Generate CSI-format index for BAM files. If input isn't cram.\n", 
                    "id": "#samtools-sort-index.cwl/csi"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 24, 
                        "prefix": "-m"
                    }, 
                    "doc": "Set minimum interval size for CSI indices to 2^INT [14]\n", 
                    "id": "#samtools-sort-index.cwl/interval"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "inputBinding": {
                        "position": 14, 
                        "prefix": "-n"
                    }, 
                    "doc": "Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates\n", 
                    "id": "#samtools-sort-index.cwl/sort_by_name"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "inputBinding": {
                        "position": 11, 
                        "prefix": "-l"
                    }, 
                    "doc": "SORT: desired compression level for the final output file, ranging from 0 (uncompressed)\nor 1 (fastest but minimal compression) to 9 (best compression but slowest to write),\nsimilarly to gzip(1)'s compression level setting.\nIf -l is not used, the default compression level will apply.\n", 
                    "id": "#samtools-sort-index.cwl/sort_compression_level"
                }, 
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 16, 
                        "valueFrom": "$(self.basename)"
                    }, 
                    "doc": "Input only in.sam|in.bam|in.cram\n", 
                    "id": "#samtools-sort-index.cwl/sort_input"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 12, 
                        "prefix": "-o", 
                        "valueFrom": "${\n    if (self == null || inputs.trigger == false){\n      return default_bam();\n    } else {\n      return self;\n    }\n}\n"
                    }, 
                    "default": null, 
                    "doc": "Write the final sorted output to FILE, rather than to standard output.\nOnly out.bam|out.cram\n", 
                    "id": "#samtools-sort-index.cwl/sort_output_filename"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "Set number of sorting and compression threads [1] (Only for sorting)\n", 
                    "id": "#samtools-sort-index.cwl/threads"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "default": true, 
                    "doc": "If true - run samtools, if false - return sort_input, previously staged into output directory\n", 
                    "id": "#samtools-sort-index.cwl/trigger"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "${\n    if (inputs.sort_output_filename == null || inputs.trigger == false){\n      return default_bam();\n    } else {\n      return inputs.sort_output_filename;\n    }\n}\n"
                    }, 
                    "secondaryFiles": "${ if (inputs.sort_input.secondaryFiles && inputs.trigger == false){ return inputs.sort_input.secondaryFiles; } else { return self.location + ext(); } }", 
                    "id": "#samtools-sort-index.cwl/bam_bai_pair"
                }
            ], 
            "baseCommand": [
                "bash", 
                "-c"
            ], 
            "arguments": [
                {
                    "valueFrom": "$(inputs.trigger)", 
                    "position": 6
                }, 
                {
                    "valueFrom": "bam", 
                    "position": 13, 
                    "prefix": "-O"
                }, 
                {
                    "valueFrom": "$(inputs.threads?inputs.threads:1)", 
                    "position": 15, 
                    "prefix": "-@"
                }, 
                {
                    "valueFrom": ";", 
                    "position": 17, 
                    "shellQuote": false
                }, 
                {
                    "valueFrom": "bash", 
                    "position": 18
                }, 
                {
                    "valueFrom": "-c", 
                    "position": 19
                }, 
                {
                    "valueFrom": "$(inputs.trigger)", 
                    "position": 21
                }, 
                {
                    "valueFrom": "$(inputs.bai?'-b':inputs.csi?'-c':[])", 
                    "position": 23
                }, 
                {
                    "valueFrom": "$(inputs.threads?inputs.threads:1)", 
                    "position": 25, 
                    "prefix": "-@"
                }, 
                {
                    "valueFrom": "${\n    if (inputs.sort_output_filename == null || inputs.trigger == false){\n      return default_bam();\n    } else {\n      return inputs.sort_output_filename;\n    }\n}\n", 
                    "position": 26
                }, 
                {
                    "valueFrom": "${\n    if (inputs.sort_output_filename == null || inputs.trigger == false){\n      return default_bam() + ext();\n    } else {\n      return inputs.sort_output_filename + ext();\n    }\n}\n", 
                    "position": 27
                }
            ], 
            "s:name": "samtools-sort-index", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "This tool is used to sort and index input BAM/SAM file by means of samtools sort/index\nInput Trigger (default: true) allows to skip all calculation and return\nall input files unchanged. To set files to be returned in case of Trigger == false,\nuse the following inputs:\n  sort_input\n", 
            "id": "#samtools-sort-index.cwl", 
            "http://schema.org/mainEntity": {
                "$import": "#samtools-metadata.yaml"
            }, 
            "http://schema.org/downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/samtools-sort-index.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:michael.kotliar@cchmc.org", 
                                            "http://schema.org/name": "Michael Kotliar", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0002-6486-3898"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "Usage: samtools sort [options...] [in.bam] Options:\n  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)\n  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]\n  -n         Sort by read name\n  -o FILE    Write final output to FILE rather than standard output\n  -T PREFIX  Write temporary files to PREFIX.nnnn.bam\n      --input-fmt-option OPT[=VAL]\n               Specify a single input file format option in the form\n               of OPTION or OPTION=VALUE\n  -O, --output-fmt FORMAT[,OPT[=VAL]]...\n               Specify output format (SAM, BAM, CRAM)\n      --output-fmt-option OPT[=VAL]\n               Specify a single output file format option in the form\n               of OPTION or OPTION=VALUE\n      --reference FILE\n               Reference sequence FASTA FILE [null]\n  -@, --threads INT\n               Number of additional threads to use [0]\n\nUsage: samtools index [-bc] [-m INT] <in.bam> [out.index] Options:\n  -b       Generate BAI-format index for BAM files [default]\n  -c       Generate CSI-format index for BAM files\n  -m INT   Set minimum interval size for CSI indices to 2^INT [14]\n  -@ INT   Sets the number of threads [none]\n"
        }, 
        {
            "class": "CommandLineTool", 
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "$import": "#envvar-global.yml"
                }
            ], 
            "hints": [
                {
                    "class": "DockerRequirement", 
                    "dockerPull": "scidap/ucsc-userapps:v325", 
                    "dockerFile": "$import: ../dockerfiles/ucsc_utils/Dockerfile\n"
                }
            ], 
            "inputs": [
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "inputBinding": {
                        "position": 4, 
                        "valueFrom": "${\n    if (self == null){\n      return inputs.input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".bigwig\";\n    } else {\n      return self;\n    }\n}\n"
                    }, 
                    "default": null, 
                    "id": "#ucsc-bedgraphtobigwig.cwl/bigWig"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "-blockSize=N - Number of items to bundle in r-tree.  Default 256\n", 
                    "inputBinding": {
                        "separate": false, 
                        "position": 1, 
                        "prefix": "-blockSize="
                    }, 
                    "id": "#ucsc-bedgraphtobigwig.cwl/blockSize"
                }, 
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 3
                    }, 
                    "id": "#ucsc-bedgraphtobigwig.cwl/genomeFile"
                }, 
                {
                    "type": "File", 
                    "inputBinding": {
                        "position": 2
                    }, 
                    "id": "#ucsc-bedgraphtobigwig.cwl/input"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "doc": "-itemsPerSlot=N - Number of data points bundled at lowest level. Default 1024\n", 
                    "inputBinding": {
                        "separate": false, 
                        "position": 1, 
                        "prefix": "-itemsPerSlot="
                    }, 
                    "id": "#ucsc-bedgraphtobigwig.cwl/itemsPerSlot"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "doc": "If set, do not use compression.", 
                    "inputBinding": {
                        "position": 1, 
                        "prefix": "-unc"
                    }, 
                    "id": "#ucsc-bedgraphtobigwig.cwl/unc"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputBinding": {
                        "glob": "${\n    if (inputs.bigWig == null){\n      return inputs.input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+\".bigwig\";\n    } else {\n      return inputs.bigWig;\n    }\n}\n"
                    }, 
                    "id": "#ucsc-bedgraphtobigwig.cwl/bigWigOut"
                }
            ], 
            "baseCommand": [
                "bedGraphToBigWig"
            ], 
            "s:name": "ucsc-bedgraphtobigwig", 
            "s:codeRepository": "https://github.com/SciDAP/workflows", 
            "s:isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "doc": "Tool is used to convert bedGraph to bigWig file\n", 
            "id": "#ucsc-bedgraphtobigwig.cwl", 
            "http://schema.org/mainEntity": {
                "class": "http://schema.org/SoftwareSourceCode", 
                "s:about": "UCSC genome browser bioinformatic utilities\n", 
                "s:targetProduct": {
                    "class": "http://schema.org/SoftwareApplication", 
                    "s:applicationCategory": "commandline tool", 
                    "http://schema.org/softwareVersion": "v325", 
                    "http://schema.org/downloadURL": "http://hgdownload.cse.ucsc.edu/admin/exe/userApps.v325.src.tgz"
                }, 
                "s:license": [
                    "https://opensource.org/licenses/GPL-3.0"
                ], 
                "id": "#ucsc-metadata.yaml", 
                "name": "#ucsc-metadata.yaml", 
                "http://schema.org/name": "UCSC userApps", 
                "http://schema.org/url": "https://genome.ucsc.edu/util.html", 
                "http://schema.org/programmingLanguage": "C++"
            }, 
            "http://schema.org/downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/tools/ucsc-bedgraphtobigwig.cwl", 
            "http://schema.org/license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "http://schema.org/creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:Andrey.Kartashov@cchmc.org", 
                                            "http://schema.org/name": "Andrey Kartashov", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0001-9102-5681"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "http://schema.org/about": "usage:\n   bedGraphToBigWig in.bedGraph chrom.sizes out.bw\nwhere in.bedGraph is a four column file in the format:\n      <chrom> <start> <end> <value>\nand chrom.sizes is a two-column file/URL: <chromosome name> <size in bases> and out.bw is the output indexed big wig file. If the assembly <db> is hosted by UCSC, chrom.sizes can be a URL like\n  http://hgdownload.cse.ucsc.edu/goldenPath/<db>/bigZips/<db>.chrom.sizes\nor you may use the script fetchChromSizes to download the chrom.sizes file. If not hosted by UCSC, a chrom.sizes file can be generated by running twoBitInfo on the assembly .2bit file. The input bedGraph file must be sorted, use the unix sort command:\n  sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph\noptions:\n   -blockSize=N - Number of items to bundle in r-tree.  Default 256\n   -itemsPerSlot=N - Number of data points bundled at lowest level. Default 1024\n   -unc - If set, do not use compression.\n"
        }, 
        {
            "class": "Workflow", 
            "requirements": [
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }
            ], 
            "doc": "creates genome coverage bigWig file from .bam file", 
            "inputs": [
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/bigWig"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/fragmentsize"
                }, 
                {
                    "type": "File", 
                    "id": "#bam-genomecov-bigwig.cwl/genomeFile"
                }, 
                {
                    "type": "File", 
                    "id": "#bam-genomecov-bigwig.cwl/input"
                }, 
                {
                    "type": [
                        "null", 
                        "double"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/mappedreads"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/pairchip"
                }, 
                {
                    "type": [
                        "null", 
                        "float"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/scale"
                }, 
                {
                    "type": [
                        "null", 
                        "boolean"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/split"
                }, 
                {
                    "type": [
                        "null", 
                        "string"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/strand"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "outputSource": "#bam-genomecov-bigwig.cwl/genomecov/genomecoverage", 
                    "id": "#bam-genomecov-bigwig.cwl/bed_file"
                }, 
                {
                    "type": "File", 
                    "outputSource": "#bam-genomecov-bigwig.cwl/bigwig/bigWigOut", 
                    "id": "#bam-genomecov-bigwig.cwl/outfile"
                }
            ], 
            "steps": [
                {
                    "run": "#ucsc-bedgraphtobigwig.cwl", 
                    "in": [
                        {
                            "source": "#bam-genomecov-bigwig.cwl/bigWig", 
                            "id": "#bam-genomecov-bigwig.cwl/bigwig/bigWig"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/genomeFile", 
                            "id": "#bam-genomecov-bigwig.cwl/bigwig/genomeFile"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/sort/sorted", 
                            "id": "#bam-genomecov-bigwig.cwl/bigwig/input"
                        }
                    ], 
                    "out": [
                        "#bam-genomecov-bigwig.cwl/bigwig/bigWigOut"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/bigwig"
                }, 
                {
                    "run": "#bedtools-genomecov.cwl", 
                    "in": [
                        {
                            "default": "-bg", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/dept"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/fragmentsize", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/fragmentsize"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/genomeFile", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/genomeFile"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/input", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/input"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/mappedreads", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/mappedreads"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/pairchip", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/pairchip"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/scale", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/scale"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/split", 
                            "valueFrom": "${\n  if (self == null){\n    return true;\n  } else {\n    return self;\n  }\n}\n", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/split"
                        }, 
                        {
                            "source": "#bam-genomecov-bigwig.cwl/strand", 
                            "id": "#bam-genomecov-bigwig.cwl/genomecov/strand"
                        }
                    ], 
                    "out": [
                        "#bam-genomecov-bigwig.cwl/genomecov/genomecoverage"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/genomecov"
                }, 
                {
                    "run": "#linux-sort.cwl", 
                    "in": [
                        {
                            "source": "#bam-genomecov-bigwig.cwl/genomecov/genomecoverage", 
                            "id": "#bam-genomecov-bigwig.cwl/sort/input"
                        }, 
                        {
                            "default": [
                                "1,1", 
                                "2,2n"
                            ], 
                            "id": "#bam-genomecov-bigwig.cwl/sort/key"
                        }
                    ], 
                    "out": [
                        "#bam-genomecov-bigwig.cwl/sort/sorted"
                    ], 
                    "id": "#bam-genomecov-bigwig.cwl/sort"
                }
            ], 
            "id": "#bam-genomecov-bigwig.cwl"
        }, 
        {
            "class": "Workflow", 
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                }, 
                {
                    "class": "ScatterFeatureRequirement"
                }, 
                {
                    "class": "StepInputExpressionRequirement"
                }, 
                {
                    "class": "InlineJavascriptRequirement"
                }, 
                {
                    "class": "MultipleInputFeatureRequirement"
                }
            ], 
            "inputs": [
                {
                    "type": "Directory", 
                    "label": "BOWTIE indices folder", 
                    "doc": "Path to BOWTIE generated indices folder", 
                    "id": "#main/bowtie_indices_folder"
                }, 
                {
                    "type": "boolean", 
                    "label": "Callpeak broad", 
                    "doc": "Set to call broad peak for MACS2", 
                    "id": "#main/broad_peak"
                }, 
                {
                    "type": "File", 
                    "label": "Chromosome length file", 
                    "format": "http://edamontology.org/format_2330", 
                    "doc": "Chromosome length file", 
                    "id": "#main/chrom_length"
                }, 
                {
                    "type": "int", 
                    "label": "Clip from 3p end", 
                    "doc": "Number of bases to clip from the 3p end", 
                    "id": "#main/clip_3p_end"
                }, 
                {
                    "type": "int", 
                    "label": "Clip from 5p end", 
                    "doc": "Number of bases to clip from the 5p end", 
                    "id": "#main/clip_5p_end"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "Control BAM file", 
                    "format": "http://edamontology.org/format_2572", 
                    "doc": "Control BAM file file for MACS2 peak calling", 
                    "id": "#main/control_file"
                }, 
                {
                    "type": "int", 
                    "label": "Expected fragment size", 
                    "doc": "Expected fragment size for MACS2", 
                    "id": "#main/exp_fragment_size"
                }, 
                {
                    "type": "File", 
                    "label": "FASTQ input file", 
                    "format": "http://edamontology.org/format_1930", 
                    "doc": "Reads data in a FASTQ format, received after single end sequencing", 
                    "id": "#main/fastq_input_file"
                }, 
                {
                    "type": "boolean", 
                    "label": "Force fragment size", 
                    "doc": "Force MACS2 to use exp_fragment_size", 
                    "id": "#main/force_fragment_size"
                }, 
                {
                    "type": "string", 
                    "label": "Effective genome size", 
                    "doc": "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9", 
                    "id": "#main/genome_size"
                }, 
                {
                    "type": "boolean", 
                    "label": "Remove duplicates", 
                    "doc": "Calls samtools rmdup to remove duplicates from sortesd BAM file", 
                    "id": "#main/remove_duplicates"
                }, 
                {
                    "type": [
                        "null", 
                        "int"
                    ], 
                    "label": "Threads", 
                    "doc": "Number of threads for those steps that support multithreading", 
                    "id": "#main/threads"
                }
            ], 
            "outputs": [
                {
                    "type": "File", 
                    "format": "http://edamontology.org/format_2572", 
                    "label": "Coordinate sorted BAM alignment file (+index BAI)", 
                    "doc": "Coordinate sorted BAM file and BAI index file", 
                    "outputSource": "#main/samtools_sort_index_after_rmdup/bam_bai_pair", 
                    "id": "#main/bambai_pair"
                }, 
                {
                    "type": "File", 
                    "format": "http://edamontology.org/format_3006", 
                    "label": "BigWig file", 
                    "doc": "Generated BigWig file", 
                    "outputSource": "#main/bam_to_bigwig/outfile", 
                    "id": "#main/bigwig"
                }, 
                {
                    "type": "File", 
                    "label": "BOWTIE alignment log", 
                    "format": "http://edamontology.org/format_2330", 
                    "doc": "BOWTIE generated alignment log", 
                    "outputSource": "#main/bowtie_aligner/output_bowtie_log", 
                    "id": "#main/bowtie_log"
                }, 
                {
                    "type": "File", 
                    "label": "FASTQ statistics", 
                    "format": "http://edamontology.org/format_2330", 
                    "doc": "fastx_quality_stats generated FASTQ file quality statistics file", 
                    "outputSource": "#main/fastx_quality_stats/statistics", 
                    "id": "#main/fastx_statistics"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "Bowtie & Samtools Rmdup combined log", 
                    "format": "http://edamontology.org/format_2330", 
                    "doc": "Processed and combined Bowtie aligner and Samtools rmdup log", 
                    "outputSource": "#main/get_stat/output", 
                    "id": "#main/get_stat_log"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "Broad peaks", 
                    "format": "http://edamontology.org/format_3614", 
                    "doc": "Contains the peak locations together with peak summit, pvalue and qvalue", 
                    "outputSource": "#main/macs2_callpeak_forced/broad_peak_file", 
                    "id": "#main/macs_broad_peaks"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "Called peaks", 
                    "format": "http://edamontology.org/format_3468", 
                    "doc": "XLS file to include information about called peaks", 
                    "outputSource": "#main/macs2_callpeak_forced/peak_xls_file", 
                    "id": "#main/macs_called_peaks"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "Gapped peak", 
                    "format": "http://edamontology.org/format_3586", 
                    "doc": "Contains both the broad region and narrow peaks", 
                    "outputSource": "#main/macs2_callpeak_forced/gapped_peak_file", 
                    "id": "#main/macs_gapped_peak"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "MACS2 log", 
                    "format": "http://edamontology.org/format_2330", 
                    "doc": "MACS2 output log", 
                    "outputSource": "#main/macs2_callpeak_forced/macs_log", 
                    "id": "#main/macs_log"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "MACS2 generated R script", 
                    "format": "http://edamontology.org/format_2330", 
                    "doc": "R script to produce a PDF image about the model based on your data", 
                    "outputSource": "#main/macs2_callpeak_forced/moder_r_file", 
                    "id": "#main/macs_moder_r"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "Narrow peaks", 
                    "format": "http://edamontology.org/format_3613", 
                    "doc": "Contains the peak locations together with peak summit, pvalue and qvalue", 
                    "outputSource": "#main/macs2_callpeak_forced/narrow_peak_file", 
                    "id": "#main/macs_narrow_peaks"
                }, 
                {
                    "type": [
                        "null", 
                        "File"
                    ], 
                    "label": "Peak summits", 
                    "format": "http://edamontology.org/format_3003", 
                    "doc": "Contains the peak summits locations for every peaks", 
                    "outputSource": "#main/macs2_callpeak_forced/peak_summits_file", 
                    "id": "#main/macs_peak_summits"
                }, 
                {
                    "type": "File", 
                    "label": "Remove duplicates log", 
                    "format": "http://edamontology.org/format_2330", 
                    "doc": "Samtools rmdup generated log", 
                    "outputSource": "#main/samtools_rmdup/rmdup_log", 
                    "id": "#main/samtools_rmdup_log"
                }
            ], 
            "steps": [
                {
                    "run": "#bam-genomecov-bigwig.cwl", 
                    "in": [
                        {
                            "source": "#main/chrom_length", 
                            "id": "#main/bam_to_bigwig/genomeFile"
                        }, 
                        {
                            "source": "#main/samtools_sort_index_after_rmdup/bam_bai_pair", 
                            "id": "#main/bam_to_bigwig/input"
                        }, 
                        {
                            "source": "#main/bamtools_stats/mappedreads", 
                            "id": "#main/bam_to_bigwig/mappedreads"
                        }
                    ], 
                    "out": [
                        "#main/bam_to_bigwig/outfile"
                    ], 
                    "id": "#main/bam_to_bigwig"
                }, 
                {
                    "run": "#bamtools-stats.cwl", 
                    "in": [
                        {
                            "source": "#main/samtools_sort_index_after_rmdup/bam_bai_pair", 
                            "id": "#main/bamtools_stats/input_files"
                        }
                    ], 
                    "out": [
                        "#main/bamtools_stats/mappedreads"
                    ], 
                    "id": "#main/bamtools_stats"
                }, 
                {
                    "run": "#bowtie.cwl", 
                    "in": [
                        {
                            "default": true, 
                            "id": "#main/bowtie_aligner/best"
                        }, 
                        {
                            "source": "#main/clip_3p_end", 
                            "id": "#main/bowtie_aligner/clip_3p_end"
                        }, 
                        {
                            "source": "#main/clip_5p_end", 
                            "id": "#main/bowtie_aligner/clip_5p_end"
                        }, 
                        {
                            "source": "#main/fastq_input_file", 
                            "id": "#main/bowtie_aligner/filelist"
                        }, 
                        {
                            "source": "#main/bowtie_indices_folder", 
                            "id": "#main/bowtie_aligner/indices_folder"
                        }, 
                        {
                            "default": 1, 
                            "id": "#main/bowtie_aligner/m"
                        }, 
                        {
                            "default": true, 
                            "id": "#main/bowtie_aligner/q"
                        }, 
                        {
                            "default": true, 
                            "id": "#main/bowtie_aligner/sam"
                        }, 
                        {
                            "default": true, 
                            "id": "#main/bowtie_aligner/strata"
                        }, 
                        {
                            "source": "#main/threads", 
                            "id": "#main/bowtie_aligner/threads"
                        }, 
                        {
                            "default": 3, 
                            "id": "#main/bowtie_aligner/v"
                        }
                    ], 
                    "out": [
                        "#main/bowtie_aligner/output", 
                        "#main/bowtie_aligner/output_bowtie_log"
                    ], 
                    "id": "#main/bowtie_aligner"
                }, 
                {
                    "run": "#fastx-quality-stats.cwl", 
                    "in": [
                        {
                            "source": "#main/fastq_input_file", 
                            "id": "#main/fastx_quality_stats/input_file"
                        }
                    ], 
                    "out": [
                        "#main/fastx_quality_stats/statistics"
                    ], 
                    "id": "#main/fastx_quality_stats"
                }, 
                {
                    "run": "#python-get-stat.cwl", 
                    "in": [
                        {
                            "source": "#main/bowtie_aligner/output_bowtie_log", 
                            "id": "#main/get_stat/bowtie_log"
                        }, 
                        {
                            "source": "#main/samtools_rmdup/rmdup_log", 
                            "id": "#main/get_stat/rmdup_log"
                        }
                    ], 
                    "out": [
                        "#main/get_stat/output"
                    ], 
                    "id": "#main/get_stat"
                }, 
                {
                    "run": "#macs2-callpeak.cwl", 
                    "in": [
                        {
                            "source": "#main/broad_peak", 
                            "id": "#main/macs2_callpeak/broad"
                        }, 
                        {
                            "default": 10000, 
                            "id": "#main/macs2_callpeak/buffer_size"
                        }, 
                        {
                            "source": [
                                "#main/force_fragment_size", 
                                "#main/exp_fragment_size"
                            ], 
                            "valueFrom": "${\n  if (!self[0]){\n    return self[1];\n  } else {\n    return null;\n  }\n}\n", 
                            "id": "#main/macs2_callpeak/bw"
                        }, 
                        {
                            "source": "#main/broad_peak", 
                            "valueFrom": "$(!self)", 
                            "id": "#main/macs2_callpeak/call_summits"
                        }, 
                        {
                            "source": "#main/control_file", 
                            "id": "#main/macs2_callpeak/control"
                        }, 
                        {
                            "source": [
                                "#main/force_fragment_size", 
                                "#main/exp_fragment_size"
                            ], 
                            "valueFrom": "${\n  if (self[0]){\n    return self[1];\n  } else {\n    return null;\n  }\n}\n", 
                            "id": "#main/macs2_callpeak/extsize"
                        }, 
                        {
                            "default": "BAM", 
                            "id": "#main/macs2_callpeak/format"
                        }, 
                        {
                            "source": "#main/genome_size", 
                            "id": "#main/macs2_callpeak/genome_size"
                        }, 
                        {
                            "default": "auto", 
                            "id": "#main/macs2_callpeak/keep_dup"
                        }, 
                        {
                            "default": "4 40", 
                            "id": "#main/macs2_callpeak/mfold"
                        }, 
                        {
                            "source": "#main/control_file", 
                            "valueFrom": "${\n  return !Boolean(self);\n}\n", 
                            "id": "#main/macs2_callpeak/nolambda"
                        }, 
                        {
                            "source": "#main/force_fragment_size", 
                            "valueFrom": "${\n  return Boolean(self);\n}\n", 
                            "id": "#main/macs2_callpeak/nomodel"
                        }, 
                        {
                            "default": 0.05, 
                            "id": "#main/macs2_callpeak/q_value"
                        }, 
                        {
                            "source": "#main/samtools_sort_index_after_rmdup/bam_bai_pair", 
                            "id": "#main/macs2_callpeak/treatment"
                        }, 
                        {
                            "default": 3, 
                            "id": "#main/macs2_callpeak/verbose"
                        }
                    ], 
                    "out": [
                        "#main/macs2_callpeak/peak_xls_file", 
                        "#main/macs2_callpeak/narrow_peak_file", 
                        "#main/macs2_callpeak/peak_summits_file", 
                        "#main/macs2_callpeak/broad_peak_file", 
                        "#main/macs2_callpeak/moder_r_file", 
                        "#main/macs2_callpeak/gapped_peak_file", 
                        "#main/macs2_callpeak/treat_pileup_bdg_file", 
                        "#main/macs2_callpeak/control_lambda_bdg_file", 
                        "#main/macs2_callpeak/macs_log"
                    ], 
                    "id": "#main/macs2_callpeak"
                }, 
                {
                    "run": "#macs2-callpeak.cwl", 
                    "in": [
                        {
                            "source": "#main/broad_peak", 
                            "id": "#main/macs2_callpeak_forced/broad"
                        }, 
                        {
                            "source": "#main/macs2_callpeak/broad_peak_file", 
                            "id": "#main/macs2_callpeak_forced/broad_peak_file_staged"
                        }, 
                        {
                            "default": 10000, 
                            "id": "#main/macs2_callpeak_forced/buffer_size"
                        }, 
                        {
                            "source": "#main/broad_peak", 
                            "valueFrom": "$(!self)", 
                            "id": "#main/macs2_callpeak_forced/call_summits"
                        }, 
                        {
                            "source": "#main/control_file", 
                            "id": "#main/macs2_callpeak_forced/control"
                        }, 
                        {
                            "source": "#main/macs2_callpeak/control_lambda_bdg_file", 
                            "id": "#main/macs2_callpeak_forced/control_lambda_bdg_file_staged"
                        }, 
                        {
                            "source": "#main/exp_fragment_size", 
                            "id": "#main/macs2_callpeak_forced/extsize"
                        }, 
                        {
                            "default": "BAM", 
                            "id": "#main/macs2_callpeak_forced/format"
                        }, 
                        {
                            "source": "#main/macs2_callpeak/gapped_peak_file", 
                            "id": "#main/macs2_callpeak_forced/gapped_peak_file_staged"
                        }, 
                        {
                            "source": "#main/genome_size", 
                            "id": "#main/macs2_callpeak_forced/genome_size"
                        }, 
                        {
                            "default": "auto", 
                            "id": "#main/macs2_callpeak_forced/keep_dup"
                        }, 
                        {
                            "default": "4 40", 
                            "id": "#main/macs2_callpeak_forced/mfold"
                        }, 
                        {
                            "source": "#main/macs2_callpeak/moder_r_file", 
                            "id": "#main/macs2_callpeak_forced/moder_r_file_staged"
                        }, 
                        {
                            "source": "#main/macs2_callpeak/narrow_peak_file", 
                            "id": "#main/macs2_callpeak_forced/narrow_peak_file_staged"
                        }, 
                        {
                            "source": "#main/control_file", 
                            "valueFrom": "${\n  return !Boolean(self);\n}\n", 
                            "id": "#main/macs2_callpeak_forced/nolambda"
                        }, 
                        {
                            "default": true, 
                            "id": "#main/macs2_callpeak_forced/nomodel"
                        }, 
                        {
                            "source": "#main/macs2_callpeak/peak_summits_file", 
                            "id": "#main/macs2_callpeak_forced/peak_summits_file_staged"
                        }, 
                        {
                            "source": "#main/macs2_callpeak/peak_xls_file", 
                            "id": "#main/macs2_callpeak_forced/peak_xls_file_staged"
                        }, 
                        {
                            "default": 0.05, 
                            "id": "#main/macs2_callpeak_forced/q_value"
                        }, 
                        {
                            "source": "#main/macs2_callpeak/treat_pileup_bdg_file", 
                            "id": "#main/macs2_callpeak_forced/treat_pileup_bdg_file_staged"
                        }, 
                        {
                            "source": "#main/samtools_sort_index_after_rmdup/bam_bai_pair", 
                            "id": "#main/macs2_callpeak_forced/treatment"
                        }, 
                        {
                            "source": [
                                "#main/force_fragment_size", 
                                "#main/macs_island_count/fragments"
                            ], 
                            "valueFrom": "${\n  return !self[0] && parseInt(self[1]) < 80;\n}\n", 
                            "id": "#main/macs2_callpeak_forced/trigger"
                        }, 
                        {
                            "default": 3, 
                            "id": "#main/macs2_callpeak_forced/verbose"
                        }
                    ], 
                    "out": [
                        "#main/macs2_callpeak_forced/peak_xls_file", 
                        "#main/macs2_callpeak_forced/narrow_peak_file", 
                        "#main/macs2_callpeak_forced/peak_summits_file", 
                        "#main/macs2_callpeak_forced/broad_peak_file", 
                        "#main/macs2_callpeak_forced/moder_r_file", 
                        "#main/macs2_callpeak_forced/gapped_peak_file", 
                        "#main/macs2_callpeak_forced/treat_pileup_bdg_file", 
                        "#main/macs2_callpeak_forced/control_lambda_bdg_file", 
                        "#main/macs2_callpeak_forced/macs_log"
                    ], 
                    "id": "#main/macs2_callpeak_forced"
                }, 
                {
                    "run": "#macs2-island-count.cwl", 
                    "in": [
                        {
                            "source": "#main/macs2_callpeak/peak_xls_file", 
                            "id": "#main/macs_island_count/input_file"
                        }
                    ], 
                    "out": [
                        "#main/macs_island_count/fragments", 
                        "#main/macs_island_count/islands"
                    ], 
                    "id": "#main/macs_island_count"
                }, 
                {
                    "run": "#samtools-rmdup.cwl", 
                    "in": [
                        {
                            "source": "#main/samtools_sort_index/bam_bai_pair", 
                            "id": "#main/samtools_rmdup/input_file"
                        }, 
                        {
                            "default": true, 
                            "id": "#main/samtools_rmdup/single_end"
                        }, 
                        {
                            "source": "#main/remove_duplicates", 
                            "id": "#main/samtools_rmdup/trigger"
                        }
                    ], 
                    "out": [
                        "#main/samtools_rmdup/rmdup_output", 
                        "#main/samtools_rmdup/rmdup_log"
                    ], 
                    "id": "#main/samtools_rmdup"
                }, 
                {
                    "run": "#samtools-sort-index.cwl", 
                    "in": [
                        {
                            "source": "#main/bowtie_aligner/output", 
                            "id": "#main/samtools_sort_index/sort_input"
                        }, 
                        {
                            "source": "#main/threads", 
                            "id": "#main/samtools_sort_index/threads"
                        }
                    ], 
                    "out": [
                        "#main/samtools_sort_index/bam_bai_pair"
                    ], 
                    "id": "#main/samtools_sort_index"
                }, 
                {
                    "run": "#samtools-sort-index.cwl", 
                    "in": [
                        {
                            "source": "#main/samtools_rmdup/rmdup_output", 
                            "id": "#main/samtools_sort_index_after_rmdup/sort_input"
                        }, 
                        {
                            "source": "#main/threads", 
                            "id": "#main/samtools_sort_index_after_rmdup/threads"
                        }, 
                        {
                            "source": "#main/remove_duplicates", 
                            "id": "#main/samtools_sort_index_after_rmdup/trigger"
                        }
                    ], 
                    "out": [
                        "#main/samtools_sort_index_after_rmdup/bam_bai_pair"
                    ], 
                    "id": "#main/samtools_sort_index_after_rmdup"
                }
            ], 
            "s:downloadUrl": "https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/scidap/run-dna-single-end.cwl", 
            "s:license": "http://www.apache.org/licenses/LICENSE-2.0", 
            "s:creator": [
                {
                    "class": "http://schema.org/Organization", 
                    "s:location": [
                        {
                            "class": "http://schema.org/PostalAddress", 
                            "s:addressLocality": "Cincinnati", 
                            "s:postalCode": "45229", 
                            "s:telephone": "+1(513)636-4200", 
                            "http://schema.org/addressCountry": "USA", 
                            "http://schema.org/addressRegion": "OH", 
                            "http://schema.org/streetAddress": "3333 Burnet Ave"
                        }
                    ], 
                    "s:department": [
                        {
                            "class": "http://schema.org/Organization", 
                            "s:department": [
                                {
                                    "class": "http://schema.org/Organization", 
                                    "s:member": [
                                        {
                                            "class": "http://schema.org/Person", 
                                            "s:email": "mailto:michael.kotliar@cchmc.org", 
                                            "http://schema.org/name": "Michael Kotliar", 
                                            "http://schema.org/sameAs": [
                                                {
                                                    "$import": "#0000-0002-6486-3898"
                                                }
                                            ]
                                        }
                                    ], 
                                    "http://schema.org/legalName": "Barski Research Lab"
                                }
                            ], 
                            "http://schema.org/legalName": "Allergy and Immunology"
                        }
                    ], 
                    "http://schema.org/legalName": "Cincinnati Children's Hospital Medical Center", 
                    "http://schema.org/logo": "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
                }
            ], 
            "doc": "Current workflow is used to run CHIP-Seq basic analysis with single-end input FASTQ file.\nIn outputs it returns coordinate sorted BAM file alongside with index BAI file, quality\nstatistics of the input FASTQ file, reads coverage in a form of BigWig file, peaks calling\ndata in a form of narrowPeak or broadPeak files.\n", 
            "id": "#main", 
            "http://schema.org/name": "run-dna-single-end", 
            "http://schema.org/codeRepository": "https://github.com/SciDAP/workflows", 
            "http://schema.org/isPartOf": {
                "class": "http://schema.org/CreativeWork", 
                "s:url": "http://commonwl.org/", 
                "http://schema.org/name": "Common Workflow Language"
            }, 
            "http://schema.org/about": "Current workflow is used to run CHIP-Seq basic analysis with single-end input FASTQ file. In outputs it returns coordinate sorted BAM file alongside with index BAI file, quality statistics of the input FASTQ file, reads coverage in a form of BigWig file, peaks calling data in a form of narrowPeak or broadPeak files. Workflow starts with running fastx_quality_stats (Step fastx_quality_stats) from FASTX-Toolkit to calculate quality statistics for input FASTQ file. At the same time (Step bowtie_aligner) Bowtie is used to align reads from input FASTQ file to reference genome. The output of this step is unsorted SAM file which is being sorted and indexed by samtools sort and samtools index (Step samtools_sort_index). Depending on workflow\u2019s input parameters indexed and sorted BAM file could be processed by samtools rmdup (Step samtools_rmdup) to remove all possible read duplicates. In a case when removing duplicates is not necessary the step returns original input BAM and BAI files without any processing. If the duplicates were removed the following step (Step samtools_sort_index_after_rmdup) reruns samtools sort and samtools index with BAM and BAI files, if not - the step returns original unchanged input files. Right after that (Step macs2_callpeak) macs2 callpeak performs peak calling. On the base of returned outputs the next step (Step macs_island_count) calculates the number and islands and estimated fragment size. If the last one is less that 80 (hardcoded in a workflow) macs2 callpeak is rerun again with forced fixed fragment size value (Step macs2_callpeak_forced). If at the very beginning it was set in workflow input parameters to run peak calling with fixed fragment size, this step is skipped and the original peak calling results are saved. The following two steps (Step bamtools_stats and bam_to_bigwig) are used to calculate coverage on the base of input BAM file and save it in BigWig format. For that purpose bamtools stats returns the number of mapped reads number which is then used as scaling factor by bedtools genomecov when it performs coverage calculation and saves it in BED format. The last one is then being sorted and converted to BigWig format by bedGraphToBigWig tool from UCSC utilities.\n"
        }
    ]
}