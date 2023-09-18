#!/usr/bin/python

import sys

with open(sys.argv[1], mode='r') as in_gff_fh:
	with open(sys.argv[2], mode='w') as out_gff_fh:

		for entry in in_gff_fh:
			if not entry.startswith("#"):
				entry_split_list = entry.split("\t")
				
				if entry_split_list[6] == "?":
					entry_split_list[6] = "."

				out_gff_fh.write("\t".join(entry_split_list))
			else:
				out_gff_fh.write(entry)
