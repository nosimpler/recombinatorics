recomb_plot <- recomb %>%
  mutate('Codon 1' = paste(AMINO_MOM, ' (', MOM_CODON, ')', sep =''),
         'Codon 2' = paste(AMINO_DAD, '(', DAD_CODON, ')', sep =''))
ggplot(recomb_plot, aes(x = `Codon 1`, y = `Codon 2`, fill = AMINO1))+
  geom_tile()+
  geom_text(aes(label = AMINO1), size = 3)+
  labs(title = 'Upper recombination')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'))+
  theme(axis.text.y = element_text(face = 'bold'))+
  scale_fill_hue(l = 50)
ggsave('C:/Users/nosim/OneDrive/Desktop/src/nosimpler.github.io/assets/images/recombination_table1.png',
       height = 20,
       width = 20
)
ggplot(recomb_plot, aes(x = `Codon 1`, y = `Codon 2`, fill = AMINO2))+
  geom_tile()+
  geom_text(aes(label = AMINO2), size = 3)+
  labs(title = 'Lower recombination')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'))+
  theme(axis.text.y = element_text(face = 'bold'))+
  scale_fill_hue(l = 50)
ggsave('C:/Users/nosim/OneDrive/Desktop/src/nosimpler.github.io/assets/images/recombination_table2.png',
       height = 20,
       width = 20
)