#df <- read.table("hic.summary", sep="\t", col.names = c("Name", "Chr1", "Pos1", "Str1", "Chr2", "Pos2", "Str2"))
df <- read.table("merged_nodups.txt", sep = " ")

df$group <- sample(rep(1:2, length.out=nrow(df)))

subdf1 <- df[df$group == 1, ]
subdf2 <- df[df$group == 2, ]

write.table(subdf1, file = "mnd-1.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(subdf2, file = "mnd-2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
