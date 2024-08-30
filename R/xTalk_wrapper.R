xTalk_wrapper <-
  function(enrich.df = NULL,
           # CrossTalk_Type = "Type_II",
           showCategory = 100,
           fc.dat   = NULL,
           doPlot = T,
           doXtalkQuant = T,
           layout = "linear",
           colorEdge = T,
           isCircular = T,
           node_label = "all",
           nPermute = 10,
           neighbourhood_th = 1,
           string_PPI_score_th = 900,
           plot_width = 10,
           plot_height = 8,
           min_cross_talk_score = 0.1,
           outdir = "./xTalk_results/") {

    set.seed(123)
    dir.create(outdir, recursive = T)
    if(!doPlot && !doXtalkQuant){
      message("Warning: doPlot and doXtalkQuant are both set as false.")
    }

    # make XtalkPlot : -----
    # fix the aesthatics: remove boundary, take care of legends, etc.
    if(doPlot){
      # spinner <- c("|", "/", "-", "\\")
      message("CrossTalk visualization (Type_I) is in progress: ")
      message("================================================\n")
      plot.obj <- xTalkPlot.Type_I.NetworkView(enrich.df = enrich.df,
                                   string_PPI_score_th = string_PPI_score_th,
                                   showCategory = showCategory,
                                   fc.dat = fc.dat,
                                   layout = layout,
                                   colorEdge = colorEdge,
                                   circular = isCircular,
                                   node_label = node_label)
      save(plot.obj, file = paste0(outdir, "xTalkPlot.Type_I.NetworkView.obj.RData"))
      plot.obj %>%
        ggsave(filename = paste0(outdir, "xTalkPlot_Type_I.pdf"),
               width = plot_width, height = plot_height, limitsize = F)
      message("Type-I Network plot [DONE]")
      xTalkPlot.Type_I.Sankey(enrich.df = enrich.df, outdir = outdir)
      message("Type-I Sankey plot [DONE]")
      message("\n================================================")

    }
    # ----

    # Make pair-wise entries for the enriched Pathways ----
    if(doXtalkQuant){
      message("CrossTalk Quantification (Type-I and Type-II) is \nin progress: ")
      message("================================================")

      # browser()
      xtalk.df <- enrich.df %>%
        full_join(enrich.df, by = character(0)) %>%
        filter(Description.x < Description.y)
      # xtalk.df <- enrich.df %>%
      #   full_join(enrich.df, by = character(0)) %>%
      #   rowwise() %>%
      #   mutate(Description1 = pmin(Description.x, Description.y),
      #          Description2 = pmax(Description.x, Description.y)) %>%
      #   filter(Description1 != Description2) %>%
      #   ungroup() %>%
      #   distinct(Description1, Description2, .keep_all = TRUE)



      res <- xtalk.df %>% pbapply(MARGIN = 1, simplify = T, FUN = function(aRow){
        p1_orgs = aRow["Overlapping.geneSymbol.x"] %>%
          str_split("/") %>%
          unlist()

        p2_orgs = aRow["Overlapping.geneSymbol.y"] %>%
          str_split("/") %>%
          unlist()

        typeI_xtalk_val = quantify_crossTalk_jaccard_index(
          neighbourhood_th = neighbourhood_th,
          string_PPI_score_th = string_PPI_score_th,
          pathway_a = p1_orgs,
          pathway_b = p2_orgs)


        typeII_xtalk_val = quantify_crossTalk_network_proximity(
          neighbourhood_th = neighbourhood_th,
          string_PPI_score_th = string_PPI_score_th,
          pathway_a = p1_orgs,
          pathway_b = p2_orgs,
          nPermute = nPermute)

        return(cbind(aRow["Description.x"],
                     aRow["Description.y"],
                     typeI_xtalk_val[1],
                     typeII_xtalk_val[1],
                     typeII_xtalk_val[2],
                     typeII_xtalk_val[3],
                     typeI_xtalk_val[2],
                     typeI_xtalk_val[3]
        ))
      }) %>% t

      colnames(res) = c("Pathway_1",
                        "Pathway_2",
                        "Type_I CrossTalk score (obs)",
                        "Type_II CrossTalk score (obs)",
                        "Type_II CrossTalk (z-score)",
                        "Type_II CrossTalk (p-value)",
                        "Type_I CrossTalk size",
                        "Type_I CrossTalk points"
                        )

      res %>% fwrite(paste0(outdir, "xtalk_Quant_result.csv"))
      message("[DONE]\n")
      message("\n===========================================================\n")

      # if(doPlot){
      #   res %>% as_tibble() %>%
      #     dplyr::select(c("Pathway_1",
      #                     "Pathway_2",
      #                     "Type_I CrossTalk score (obs)")) %>%
      #     dplyr::rename("weight" = "Type_I CrossTalk score (obs)") %>%
      #     dplyr::filter(weight != 0) %>%
      #     dplyr::filter(weight >= min_cross_talk_score) %>%
      #     xTalkPlot.plain(df = ., layout = layout) %>%
      #     ggsave(filename = paste0(outdir, "xTalkPlot.pdf"),
      #            width = plot_width, height = plot_height, limitsize = F)
      # }

      if(doPlot){
        message("CrossTalk visualization (Type-II) is in progress: ")
        message("\n================================================\n")
        typeII.plot.df <- res %>% as_tibble() %>%
          select(c("Pathway_1",
                          "Pathway_2",
                          "Type_II CrossTalk score (obs)")) %>%
          rename("weight" = "Type_II CrossTalk score (obs)") %>%
          mutate(weight = weight %>% as.numeric) %>%
          # dplyr::filter(weight != 0) %>%
          filter(weight >= min_cross_talk_score) %>%
          mutate(across(c("weight"), as.numeric))

        plot.obj <- typeII.plot.df %>%
          xTalkPlot.Type_II.NetworkView(df = .,
                                        isCircular = isCircular,
                                        layout = layout)
        save(plot.obj, file = paste0(outdir, "xTalkPlot.Type_II.NetworkView.obj.RData"))
        plot.obj %>%
          ggsave(filename = paste0(outdir, "xTalkPlot_Type_II.pdf"),
                 width = plot_width, height = plot_height, limitsize = F)
        message("Type-II Network plot [DONE]")

        # here put heatmap code
        # motivation: it gives clustering view of this type-II cross-talks
        # png(file = paste0(outdir, "xtalkHeatMap.png"),
        #     width = plot_width, height = plot_height,
        #     res = 1200,
        #     units = "in")
        # typeII.plot.df %>%
        #   fwrite(paste0(outdir, "xtalk_typeII.plot.df.csv"))
        pdf(file = paste0(outdir, "xtalkHeatMap.pdf"),
                width = plot_width, height = plot_height)
        # par(
        #   mar      = c(5, 5, 2, 2),
        #   xaxs     = "i",
        #   yaxs     = "i",
        #   cex.axis = 2,
        #   cex.lab  = 2
        # )
        typeII.plot.df %>%
          xTalkPlot.Type_II.heatmap()
        dev.off()
        # typeII.plot.df %>%
        #     xTalkPlot.Type_II.heatmap() %>%
        #   save_pdf(filename = paste0(outdir, "xtalkHeatMap.png"),
        #            width = plot_width,
        #            height = plot_height)

        message("Type-II Network plot [DONE]")

        message("\n================================================\n")
        message(paste0("All the results are saved in the directory named: ", outdir))
        message("\n================================================\n")
      }
    }
    # ----
  }
# enrich.df = readRDS("C:\\Users\\Hp\\OneDrive - Imam university\\Documents\\Az projects\\COVID_CDR Projects\\COVID_CDR_Cyclosporine_Selinexor_comb_therapy\\data\\DB00091_KEGG_2021_Human_all.rds")
# enrich.df <- enrich.df %>%
#   dplyr::filter(p.adjust < 0.001)
# enrich.df %>%
#   xTalk_wrapper(doPlot = T,
#                 doXtalkQuant = T,
#                 nPermute = 2,
#                 min_cross_talk_score = 1.0,
#                 plot_width = 10,
#                 plot_height = 10)
