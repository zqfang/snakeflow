
library(infercnv)
library(Seurat)

WORKSPACE = "/data/bases/fangzq/20200505_HCC/"
DATASETS = "hcc.hashtag.demux.integrated.results.Robj"
FIGURES = paste0(WORKSPACE, "figures/")
# dir.create(paste0(WORKSPACE, "CNVs"), recursive = T)

setwd(WORKSPACE)

load(paste0(WORKSPACE, DATASETS))
DefaultAssay(hcc.combined) <- "RNA"
counts_matrix = GetAssayData(hcc.combined, slot="counts")


noinfercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=paste0(WORKSPACE, "CNVs/cellAnnotations.txt"),
                                    delim="\t",
                                    gene_order_file=paste0(WORKSPACE,"CNVs/gene_ordering_file.txt"),
                                    ref_group_names=c("N")# normal cells 
                                    ) 
# 这一步的一个关键参数是ref_group_name, 用于设置参考组。
# 假如你并不知道哪个组是正常，哪个组不正常，那么设置为ref_group_name=NULL, 
# 那么inferCNV会以全局平均值作为基线，这适用于有足够细胞存在差异的情况。
# 此外，这里的raw_count_matrix是排除了低质量细胞的count矩阵。

infercnv_obj = infercnv::run(noinfercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=paste0(WORKSPACE, "CNVs"), 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             num_threads=20,
                             HMM=TRUE)
# 关键参数是cutoff, 用于选择哪些基因会被用于分析（在所有细胞的平均表达量需要大于某个阈值）。
# 这个需要根据具体的测序深度来算，官方教程建议10X设置为0.1，smart-seq设置为1。
# 你可以先评估下不同阈值下的保留基因数，决定具体值。
# cluster_by_groups用于声明是否根据细胞注释文件的分组对肿瘤细胞进行分群。

save(infercnv_obj, file=paste0(WORKSPACE, "CNVs/hcc.infercnv.final.Robj"))

# 提取信息
# inferCNV会输出一个" map_metadata_from_infercnv .txt"文件用于记录每个细胞的元信息，
# 所有信息都可以从该文件中进行提取。或者使用infercnv::add_to_seurat将信息直接增加到原来的seurat对象中。
hcc.infercnv = infercnv::add_to_seurat(hcc.combined, paste0(WORKSPACE, "CNVs"))
save(hcc.infercnv, file=paste0(WORKSPACE, "hcc.hashtag.demux.integrated.infercnv.results.Robj"))

# Apply a median filtering to the expression matrix within each tumor bounds
# infercnv_obj <- infercnv::apply_median_filtering(infercnv_obj)


# plot result object
plot_cnv(infercnv_obj, out_dir = paste0(WORKSPACE, "figures/"), title = "inferCNV",
         obs_title = "Observations (Cells)", ref_title = "References (Cells)",
         cluster_by_groups = TRUE, cluster_references = TRUE,
         #plot_chr_scale = FALSE, chr_lengths = NULL, 
         k_obs_groups = 3,
         contig_cex = 1, x.center = mean(infercnv_obj@expr.data),
         x.range = "auto", hclust_method = "ward.D",
         #custom_color_pal = NULL, # Has to be in the shape color.palette(c("darkblue","white", "darkred"), c(2, 2))
         color_safe_pal = FALSE,
         output_filename = "infercnv",
         output_format = "png", png_res = 300,
         dynamic_resize = 0, ref_contig = NULL, write_expr_matrix = FALSE,
         useRaster = TRUE)

# Takes an infercnv object and subdivides it into one object per group of cells to allow plotting of each
# group on a seperate plot.
plot_per_group(infercnv_obj, on_references = TRUE,
               on_observations = TRUE, sample = FALSE, n_cells = 1000,
               every_n = NULL, above_m = 1000,
               base_filename = "infercnv_per_group", output_format = "png",
               write_expr_matrix = TRUE, save_objects = FALSE, png_res = 300,
               dynamic_resize = 0, out_dir)