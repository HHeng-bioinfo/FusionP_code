# ============================================================
# Sec 1: Library and Load
# ============================================================

# Sec 1.1: Libraries
library(circlize)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(grid)
library(IRanges)
library(paletteer)

set.seed(123)




##### __Sec 1.2 load #####
{
  # OMAP result
  pls_dt <- fread("OMAP_pls.tsv")
  pls_raw <- fread("OMAP_pls_raw.tsv")
  
  # blast result
  blast <- fread("blast.out")
  
  # mob result
  mob <- fread("mob_type.tsv")
  
  # bakta result
  anno <- fread("bakta_anno.tsv")
  
  # sorted data
  bact_meta <- fread("bact_meta.tsv")
  
  # meta data
  meta <- fread("PDG000000012.1904.metadata.tsv")
  
  # ref clu
  pls_clu <- fread("plasmids.repr.clu")
  
  # resistenc meta
  res_type <- fread("CARD_AMR_clustered.csv") # from kleborate
  
  
  ##
  out_dir <- "./"
  dir.create(out_dir, recursive = T)
  # 
}



# ============================================================
# Sec 2: Functions
# ============================================================

# Reduce length
red_len <- function(start, end) {
  rd <- reduce(IRanges(start, end))
  return(sum(rd@width))
}

# Calculate overlap
calculate_overlap <- function(start1, end1, start2, end2) {
  overlap <- pmax(0, pmin(end1, end2) - pmax(start1, start2))
  overlap_percentage <- overlap / pmax(end1 - start1, end2 - start2)
  return(overlap_percentage)
}



# ============================================================
# Sec 3: Meta modify and filter
# ============================================================

# Sec 3.1: Meta 
all(sort(unique(blast[!grepl("_PC_", filename), filename])) == sort(unique(blast[grepl("__", filename), filename])))
setdiff(unique(blast[!grepl("_PC_", filename), filename]), unique(blast[grepl("__", filename), filename]))
setdiff(unique(blast[grepl("__", filename), filename]), unique(blast[!grepl("_PC_", filename), filename]))
db_fn <- sort(unique(blast[!grepl("_PC_|_VC_", filename), filename]))

# CR genes
{
  CR_gene <- unique(gsub("\\.v[1-9]$", "", res_type[bla_class == "Bla_Carb", queryID]), res_type[bla_class == "Bla_Carb", queryID])
}



# MOB normalization
{
  mob <- mob[, .(sample_id, rep_type = `rep_type(s)`, predicted_mobility, size)]
  mob[rep_type == "-", rep_type := ""]
}


# pls_dt processing
print(paste0("Start: pls_dt "))
pls_dt <- pls_dt[strain %in% bact_meta$strain]
pls_dt <- pls_dt[!grepl("VC", cluster)]
qstrain <- setdiff(pls_dt[grepl("-", strain), unique(strain)], "pls__pKP00-9T")

# pls_clu reshape
pls_clu <- pls_clu[, .(strain = V1, pls_id = V1, cluster = V3, scovs = 1, slen = V4)]
pls_clu <- rbind(pls_clu, data.table(strain = "pls__pKP00-9T", pls_id = "pls__pKP00-9T", cluster = "PC_499", scovs = 1, slen = 343950))


# pls_dt processing2
pls_dt[, sample_id := ifelse(strain %in% qstrain, paste0(strain, "_", cluster), strain)]
pls_dt <- merge(pls_dt, mob, by = c("sample_id"))
pls_dt[, sample_id := NULL]
pls_dt <- merge(pls_dt, bact_meta[, .(strain, ST, KL, Year, Area)], by = "strain", all.x = T)

#
print(paste0("nrow: ", nrow(pls_dt)))

#
cols = setdiff(colnames(pls_dt), c("scovs", "slen", "Year", "size"))
pls_dt[, (cols) := lapply(.SD, function(x) ifelse(is.na(x), "", x)), .SDcols = cols]
pls_dt[, Year := as.numeric(Year)]


#
pls_dt[, filename := paste0(strain, "_", cluster)]



# Sec 3.2: BLAST filtering and annotation
#
blast <- blast[length >= 1000][pident >= 85]
print(paste0("blast filter nrow: ", nrow(blast)))
blast[filename %in% db_fn, strain := filename]
blast[filename %in% db_fn, cluster := pls_clu[, .(filename = strain, cluster)][.SD, on = "filename", x.cluster]]

#
blast[!filename %in% db_fn, strain := gsub("_PC_[0-9]*$", "", filename)]
blast[!filename %in% db_fn, cluster := sub(".*_(PC_\\d+)$", "\\1", filename)]

#
blast[, filename := paste0(strain, "_", cluster)]

#
print(paste0("Start: with blast pls clu"))
blast <- merge(blast, 
               rbind(unique(pls_raw[, .(sseqid = qseqid, s_strain = strain, s_cluster = cluster)]), 
                     pls_clu[, .(sseqid = strain, s_strain = strain, s_cluster = cluster)]), 
               allow.cartesian=TRUE)

#
print(paste0("blast merge nrow: ", nrow(blast)))

#
blast[, filename_s := paste0(s_strain, "_", s_cluster)]
blast <- blast[filename_s %in% unique(blast$filename)] ### !!!

#
print(paste0("Start: with blast pls rep_type predicted_mobility"))
blast <- merge(blast, pls_dt[, .(filename, q_size = size, rep_type, predicted_mobility)], by = c("filename"))
blast <- merge(blast, pls_dt[, .(filename_s = filename, s_size = size, rep_type_s = rep_type, predicted_mobility_s = predicted_mobility)], by = c("filename_s"))

#
print(paste0("blast merge nrow: ", nrow(blast)))

#
blast <- blast[q_size < s_size]

#
print(paste0("blast filter q/slen nrow: ", nrow(blast)))

#
blast_dt <- copy(blast)
blast_dt[send <= sstart, c("sstart", "send") := .(send, sstart)]

# 
fwrite(blast_dt, paste0(out_dir, "/", "blast_dt.tsv"), sep = '\t', eol = '\n')



# Sec 3.3: BLAST map
## blast
{
  #
  print(paste0("Start: with blast map"))
  
  #
  blast_map <- blast_dt[, red_len(sstart, send), by = .(strain, s_strain, cluster, s_cluster)]
  blast_map <- merge(blast_map, unique(blast_dt[, .(strain, cluster, q_size)]), by = c("strain", "cluster"))
  blast_map <- merge(blast_map, unique(blast_dt[, .(s_strain, s_cluster, s_size)]), by = c("s_strain", "s_cluster"))
  blast_map <- blast_map[, .(strain, s_strain, cluster, s_cluster, s_size, q_size, qcovs_pls = V1/q_size, scovs_pls = V1/s_size)]
  
  #
  fwrite(blast_map, paste0(out_dir, "/", "blast_map.tsv"), sep = '\t', eol = '\n')
}
blast_map <- fread(paste0(out_dir, "/", "blast_map.tsv"))

#
# oneseq_pls <- pls_raw[, uniqueN(qseqid), .(strain, cluster)][V1 == 1]
complete_genome <- meta[asm_level %in% c("Complete Genome", "Chromosome"), gsub("\\.", "_", gsub("GCA_", "", asm_acc))]
complete_genome <- c(complete_genome, db_fn)
rm(meta)

#
oneseq_pls <- pls_raw[, uniqueN(qseqid), .(strain, cluster)][V1 == 1]
oneseq_pls <- rbind(oneseq_pls, pls_clu[, .(strain, cluster, V1 = 1)])
gc()


# Sec 3.4: Anno processing
## colname anno
colnames(anno) <- c("#Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", 
                    "Gene", "Product", "DbXrefs")
anno_tsv <- anno[, .(locus = `Locus Tag`, name = Gene, start = Start, end = Stop, strand = Strand, product = Product)]

#
anno_tsv[, filename := gsub("_[0-9]*$", "", locus)]

#
anno_tsv[filename %in% db_fn, strain := filename]
anno_tsv[filename %in% db_fn, cluster := pls_clu[, .(filename = strain, cluster)][.SD, on = "filename", x.cluster]]

#
anno_tsv[!filename %in% db_fn, strain := gsub("_PC_[0-9]*$", "", filename)]
anno_tsv[!filename %in% db_fn, cluster := sub(".*_(PC_\\d+)$", "\\1", filename)]

#
anno_tsv[, type := "others"]

#
anno_tsv[, name := gsub("^bla", "", name)]
anno_tsv[grepl("mph", name), name := gsub("\\(|\\)", "", name)]
anno_tsv[name %in% unique(gsub("\\.v[1-9]$", "", res_type[, queryID]), res_type[, queryID]), type := "resistance"]

#
anno_tsv[product == "origin of replication", name := "oriC"]
anno_tsv[product == "origin of transfer", name := "oriT"]
anno_tsv[grepl("^rep|^mob|^ori|^tra|^trb", name, ignore.case = T), type := "plasmid_gene"]

#
anno_tsv[grepl("^iuc|^iro|^peg|^rmpA", name, ignore.case = T), type := "virulence"]


#
anno_tsv[grepl("^tnp", name, ignore.case = T), type := "IS"]
anno_tsv[grepl("^tnp", name), sub(".*\\b(\\w+) family\\b.*", "\\1", product)]
anno_tsv[grepl("^tnp", name), name := sub(".*\\b(\\w+) family\\b.*", "\\1", product)]

#
gc()


# ============================================================
# Sec 4: Data overview
# ============================================================
# Sec 4.1: BLAST map complete
#
blast_map[qcovs_pls >= 0.8]

blast_map_f <- merge(blast_map[qcovs_pls >= 0.8][s_strain %in% complete_genome], 
                     oneseq_pls[, .(s_strain = strain, s_cluster = cluster)])

fwrite(blast_map_f, paste0(out_dir, "/", "blast_map_f.tsv"), sep = '\t', eol = '\n')


##
#
blast_map_hv <- merge(blast_map_f, pls_dt[hv != "", .(s_strain = strain, s_cluster = cluster)], by = c("s_strain", "s_cluster"))
blast_map_hv <- blast_map_hv[s_cluster %in% c("PC_499")]

#
blast_map_hv[q_size >= 20000][!cluster %in% c("PC_3920", "PC_6321", "PC_6140")][!cluster %in% c("PC_499", "PC_341"), .N, cluster] ###
blast_map_hv[q_size >= 20000][, uniqueN(s_strain), cluster]

##
#
blast_map_hv <- blast_map_hv[q_size >= 20000][!cluster %in% c("PC_3920", "PC_6321", "PC_6140", "PC_5790")]


# ============================================================
# Sec 5: Plot
# ============================================================
# Sec 5.1: Plot functions
# Draw gene arrows
draw_gene_arrow <- function(start, end, strand, color) {
  if (strand == "+") {
    circos.arrow(
      start, end, 
      width = 0.8,
      arrow.head.length = 0.25,
      arrow.head.width = 0.1,
      col = color
    )
  } else if (strand == "-"){
    circos.arrow(
      end, start, 
      width = 0.8,
      arrow.head.length = 0.25,
      arrow.head.width = 0.1,
      col = color
    )
  } else {
    circos.arrow(
      end, start, 
      width = 0.8,
      arrow.head.length = 0,
      arrow.head.width = 0,
      col = color
    )
  }
}


# Add track labels in gaps
add_track_label <- function(label) {
  circos.text(0, 0, label, facing = "downward", niceFacing = T, adj = c(1.1, -1), cex = 0.8)
}


color_mapping <- colorRamp2(1994:2024, c(c(rev(paletteer_c("grDevices::Viridis", 60)))[seq(1, 20, by = 1)], 
                                         c(rev(paletteer_c("grDevices::Viridis", 30)))[seq(10, 30, by = 2)]))






# Sec 5.2: Plot hv loop
#
{
  ## #####!!!
  blast_map_in <- blast_map_hv[strain %in% blast_map_hv[cluster %in% c("PC_499"), strain]][strain %in% qstrain]
  
  ##
  out_evi <- data.table(strain = character(0), sumlen = numeric(0), rep_mb = character(0), clu = character(0), Year = numeric(0), ST = character(0), KL = character(0), s_strain = character(0), s_clu = character(0))
  blast_dt_list <- list()
  
  
  
  
  for (s_strain_in in unique(blast_map_in[!cluster %in% c("PC_499"), s_strain])) {   ### exclude 499-only mapping s_strain
    ##
    print(paste0("Start: with Ref ", s_strain_in))
    
    for (s_cluster_in in unique(blast_map_in[s_strain == s_strain_in, s_cluster])) {
      ##
      print(paste0("And Cluster: ", s_cluster_in))
      
      ##
      #
      blast_map_s_in <- blast_map_in[s_strain == s_strain_in][s_cluster == s_cluster_in]
      
      # #####!!!
      blast_map_s_in <- blast_map_s_in[strain %in% blast_map_s_in[cluster %in% c("PC_499"), strain]][strain %in% qstrain]
      
      ##
      if (nrow(blast_map_s_in) == 0) {
        next
      }
      
      ##
      #
      n_querystrain <- max(15, (blast_map_s_in[, uniqueN(strain)] + 1))
      
      #
      print(paste0("N query: ", n_querystrain))
      
      
      ##
      #
      blast_dt_in <- merge(blast_dt, blast_map_s_in[, .(s_strain, s_cluster, strain, cluster)], by = c("strain", "cluster", "s_strain", "s_cluster"))
      blast_dt_in <- merge(blast_dt_in, pls_dt[, .(strain, cluster, ST, KL, Year, Area)], by = c("strain", "cluster"))
      
      #
      if (nrow(blast_dt_in) == 0) {
        next
      }
      
      ##
      anno_tsv_in <- anno_tsv[strain == s_strain_in][cluster == s_cluster_in]
      
      
      ##
      title_message <- unique(blast_dt_in[!cluster %in% c("PC_499"), .(cluster, s_cluster, rep_type)][order(cluster, rep_type)])[
        , .(rep = paste0(rep_type, collapse = ' & '), clu = paste0(cluster, collapse = ' & '), sclu = paste0(unique(s_cluster), collapse = ' & '))]
      

      
      ## #####!!!
      pdf(paste0(out_dir, "/", "hv_", 
                 gsub(" ", "", title_message[1, sclu]), "-",  
                 gsub(" ", "", title_message[1, clu]), "_", s_strain_in, ".pdf"), width = 20, height = 20)
      {
        
        ## Set up circular genome - FIXED INITIALIZATION
        #
        circos.clear()
        circos.initialize(factors = "reference", xlim = c(0, blast_dt_in[, s_size][1]))
        
        
        ## Track 1: Gene annotation of reference plasmid 
        {
          # label
          #
          if (nrow(anno_tsv_in) > 0) {
            circos.genomicLabels(anno_tsv_in[type != "others"][name != ""][, .(chr = "reference", start, end, name)],
                                 labels.column = 4, side = "outside", cex = 0.7, connection_height = mm_h(2), 
                                 font = anno_tsv_in[type != "others"][name != ""][, ifelse(type == "IS", 1, 3)])
          } else {
            print(paste0("Annotation missing: ", s_strain_in, "_", s_cluster_in))
          }
          
          
          # anno 
          circos.track(
            ylim = c(0, 1),
            track.height = 1/n_querystrain, bg.border = NA, ### !!!
            track.margin=c(0,0), cell.padding = c(0, 0, 0, 0),
            panel.fun = function(x, y) {
              
              #
              circos.genomicAxis(h = "bottom", direction = "inside")
              
              # Draw genes as arrows
              
              if (nrow(anno_tsv_in) > 0) {
                for (i in 1:nrow(anno_tsv_in)) {
                  gene <- anno_tsv_in[i]
                  
                  # Select color based on gene type
                  if (gene$type == "resistance") {
                    gene_color = "#00468B"
                  } else if (gene$type == "virulence") {
                    gene_color = "#ED0000"
                  } else if (gene$type == "plasmid_gene") {
                    gene_color = "#78B056"
                  } else if (gene$type == "IS") {
                    gene_color = "#B254A5"
                  } else {
                    gene_color = "grey70"
                  }
                  
                  # Draw gene arrow
                  draw_gene_arrow(
                    gene$start, gene$end, gene$strand,
                    color = gene_color
                  )
                  
                }
              }
              
              
              
            }
          )
        }
        
        
        ## Track others: label prepare
        ##
        {
          blast_fus_in <- merge(unique(blast_dt_in[, .(strain, cluster, qlen)])[, .(sumlen = sum(qlen)), strain], 
                                unique(blast_dt_in[, .(strain, cluster, rep_type)][order(strain, cluster)])[
                                  , .(rep_mb = paste0(rep_type, collapse = '&'), clu = paste0(cluster, collapse = '&'))
                                  , by = .(strain)], by = ("strain"))[order(clu, -sumlen)]
          
          blast_fus_in <- merge(blast_fus_in, unique(pls_dt[, .(strain, Year, ST, KL)]))
          
          evidence_out <- rbind(blast_fus_in[(grepl("&", rep_mb))], 
                                merge(pls_dt[predicted_mobility != "non-mobilizable"], 
                                      blast_map_s_in[scovs_pls >= 0.8, .(strain, cluster)], by = c("strain", "cluster"))[
                                        , .(strain, sumlen = slen, rep_mb = rep_type, clu = cluster, Year, ST, KL)]) 
          
          evidence_out[, s_strain := s_strain_in]
          evidence_out[, s_clu := s_cluster_in]
          
          out_evi <- rbind(out_evi, evidence_out)
          # rm(evidence_out)
          
          
          ##
          blast_dt_list[[length(blast_dt_list) + 1]] <- blast_dt_in
          
          

          ##
          #
          fusion_message <- evidence_out[, paste(Year, clu, rep_mb, ST, KL, sep = '-')]
          rm(evidence_out)
          
          fusion_message <- sort(fusion_message)
          
          
          
          if (uniqueN(fusion_message) == 0) {
            fusion_message <-  "No evidence"
          }
        }
        
        
        ## Track others: mapping
        ##
        margin_i <- 0
        for (q_strain_in in blast_fus_in[order(-Year, -sumlen)]$strain) {
          print(paste0("Map: ", q_strain_in))
          
          blast_circle_in <- blast_dt_in[strain == q_strain_in]
          
          top_margin <- ifelse(margin_i == 0, 0.04, 0)
          margin_i = 1
          
          circos.track(
            ylim = c(0, 1),
            track.height = 0.7/n_querystrain,
            track.margin=c(0, top_margin), cell.padding = c(0, 0, 0, 0),
            bg.border = '#E5E5E5',
            panel.fun = function(x, y) {
              # Draw BLAST alignments
              for (i in 1:nrow(blast_circle_in)) {
                # print(i)
                alignment <- blast_circle_in[i]
                circos.rect(
                  alignment$sstart, 0.1,
                  alignment$send, 0.9,
                  col = color_mapping(alignment$Year),
                  border = NA
                )
              }
              
            }
          )
        }
        
        
        
        ## Track others: label plot
        {
          ## gene type legend
          legend(
            "bottomleft",
            title = paste0("Ref: ", pls_dt[strain == s_strain_in][cluster == s_cluster_in, 
                                                                  paste(strain, rep_type, sep = ', ')]),
            legend = c("Resistance", "Virulence", "Plasmid gene", "IS", "Other"),
            fill = c("#00468B", "#ED0000", "#78B056", "#B254A5", "grey70"),
            border = NA,
            cex = 0.8, ncol = 2
          )
          
          
          ##
          legend(
            "bottomright",
            title = paste0("Fusion evidence"),
            legend = fusion_message, 
            border = NA,
            cex = 0.8
          )
          
          
          ##
          title(main = paste0("Fusion: ", title_message[1, sclu], " and ", title_message[1, clu], " with ", title_message[1, rep]))
          
          
        }
        
      }
      dev.off()
      
      gc()
      
    }
    
    
  }
  
}



# ============================================================
# Sec 6: Output evidence
# ============================================================
#
fwrite(out_evi, paste0(out_dir, "/", "evidence.txt"), sep = '\t', eol = '\n', col.names = F)

#
if (file.exists(paste0(out_dir, "/", "blast_dt_list.txt"))) file.remove(paste0(out_dir, "/", "blast_dt_list.txt"))
fwrite(rbindlist(blast_dt_list), paste0(out_dir, "/", "blast_dt_list.txt"), sep = '\t', eol = '\n', col.names = F)


