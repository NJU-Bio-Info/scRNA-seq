#!/bin/bash

cellranger count --id <ID> \
	--transcriptome <PATH> \
	--fastqs <PATH> \
	--sample <PREFIX> \
	--lanes <NUMS> \
	--expect-cells <NUM> \
	--chemistry <CHEM> \
	--localcores <NUM> \
	--localmem <NUM>

###options
--transcriptome		the path to genome reference.
--fastqs		the path [dir] containing fastq files.
--sample		the prefix of fastq files to read in.
--lanes			the lanes to read in.
--expect-cells		the number of cells you expect to have.
--chemistry		specify the chemistry type.by default the assay configuration is detected automatically, 
			which is the recommened mode. You usually will not need to specify a chemistry. 
			Options are: 
				'auto' for autodetection
				'threeprime' for Single Cell 3'
				'fiveprime' for  Single Cell 5'
				'SC3Pv1' or 'SC3Pv2' or'SC3Pv3' for Single Cell 3' v1/v2/v3
				'SC3Pv3LT' for Single Cell 3' v3 LT
				'SC5P-PE' or 'SC5P-R2' for Single Cell 5', paired-end/R2-only
				'SC-FB' for Single Cell Antibody-only 3' v2 or 5' [default: auto]
--localcores
--localmem


###Question:
What to do if I sequenced more bases than required on the R1, R2, or the index read for 3' single-cell data?

Answer: For single indexed libraries, if you sequenced more bases than necessary on any of your reads (R1, R2, or index), you have three options.
1. You can leave the sequences "as-is." It may be helpful to specify the --chemistry option (e.g. --chemistry=SC3Pv3).
For Single Cell 3' v3.1, v3, and v2 chemistries, extra cycles on the R1 (cell-barcode, UMI) read will be ignored. By default, all of R2 (the RNA read) is used.

2. You or your sequencing core can demultiplex your BCL data using the --use-bases-mask option in cellranger mkfastq or bcl2fastq to mask the extra bases from ever appearing in your FASTQ files, which will save some disk space.

3. Alternatively, if you are unable to demultiplex again, yet you want to save as much disk space as possible, you could use a third-party tool to trim your FASTQs directly.
