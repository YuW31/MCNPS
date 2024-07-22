ddd = read.csv("n30_1.csv")
ddd = ddd[,-1]

p =   ggplot(data=ddd, aes(x=`Test.Length`, y=PAR, group=Method)) +
  # ggtitle(paste0("N=",Matrix.con[i,1]))+
  geom_line(aes(color=Method))+
  geom_point(size = 0.1) +
  scale_x_continuous(name = "Test Length", breaks = c(3,5,10,15)) +
  scale_y_continuous(name = "PAR", breaks = seq(0, 1, 0.1)) + 
  scale_color_manual(values=c('red', 'green', 'black','blue'))

p2 = p + facet_wrap(vars(Attribute,Options))

ggsave(
  "April_xxx_302.png",
  plot = p2,
  width = 5,
  height = 3,
  bg = NULL
)
