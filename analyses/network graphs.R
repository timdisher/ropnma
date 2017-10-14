library(ggnetwork)
library(sna)


pub_netgraph = function(chars, data, nodecolour = "mediumorchid2", layout = "circle"){

chars_net = chars$int_char %>% select(Treatment, `Number of Patients`) %>% rename(sample = `Number of Patients`)

df = chars$direct %>% rename(weight = `# Studies`) %>% select(Treatment, weight) %>% separate(Treatment,c("from","to"), sep = " vs ")
df = graph.data.frame(df)
V(df)$size = as.numeric(chars_net[match(V(df)$name,chars_net$Treatment),2][[1]])

legend = data.frame(treat = V(df)$name, letter = LETTERS[as.factor(V(df)$name)]) %>% mutate(legend = paste(letter,": ",treat,sep = "")) %>% select(legend) %>%
  arrange(legend)

network = ggnetwork(df, layout = layout)

ggplot(network, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(size = weight), color = "grey50", alpha = 0.8) +
  geom_nodes(aes(size = size,colour = vertex.names)) + 
  theme_blank() + 
  geom_edgelabel_repel(aes(label = weight), color = "grey25") +
  geom_nodetext(aes(label = LETTERS[vertex.names]))+
  scale_size_area("sample", max_size = 32,guide = FALSE) + scale_colour_manual(name = "Treatment", values = rep(nodecolour,length(network$x)),
                                                        labels = legend[[1]]) + theme(legend.position = "bottom")
  

}

pub_netgraph(chars = os_reac$chars, data = os_reac$data)
