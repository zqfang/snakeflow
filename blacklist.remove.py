
import os, sys


with open("temp/blacklist.txt") as black:
    blacklist = [ bla.strip().split("/")[-1] for bla in black]

#diff_Ctrl36_vs_Ctrl24/diff_Ctrl36_vs_Ctrl24_results.annotated.xls
bl = [ b.rpartition("_")[0] for b in blacklist ]
RMATS_ = [b.replace("diff_", "alternative_splicing/rMATS.") for b in bl ]
RMATS_sig = [b+"_sig" for b in RMATS_ ]
DIFF_ = ["differential_expression/"+b for b in bl ]
GO_GSEA = [ b.replace("diff", "GO/GSEA") for b in bl ]
GO_Enrichr= [ b.replace("diff", "GO/Enrichr") for b in bl ]

bl2rm = DIFF_ + GO_GSEA + GO_Enrichr + RMATS_ + RMATS_sig

for bl in bl2rm:
    os.system(" rm -r %s"%bl )

os.system("echo 'done' ")




