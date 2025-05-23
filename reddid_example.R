pacman::p_load(conflicted, tidyverse, wrappedtools, here, readxl,
               ggrepel, ggh4x)
importdata <- read_excel(here("Data/reddit_AskEurope_flair_counts.xlsx"))
colnames(importdata)[1] <- "Topic"
rawdata <- pivot_longer(importdata,
                        cols = -Topic,
                        names_to = "Date",
                        values_to = "Posts") |>
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) |>
bind_rows(
  rawdata |>
    group_by(Date) |>
    summarize(Posts=sum(Posts)) |>
    mutate(Topic="All topics")
) |>
  mutate(Combined=case_when(Topic=="All topics"~"Combined",
                            .default="Topics"))
ggplot(rawdata,aes(Date,Posts,color=Topic))+
  geom_line()+
  geom_text_repel(data=rawdata |> filter(Date=="2024-12-30"),
             aes(label=Topic),
             max.overlaps = 20)+
  scale_x_date(#date_breaks="14 days",
    # limits = c(as.Date("2024-09-23"),
    #            as.Date("2024-12-30")),
    breaks=unique(rawdata$Date),
    minor_breaks=NULL,
    guide=guide_axis(n.dodge=2),
    expand = expansion(add=c(5,10)))+
  facet_grid(rows=vars(Combined),
             scales="free_y")+
  force_panelsizes(rows = c(1, 3)) +
  theme(legend.position = "none")
