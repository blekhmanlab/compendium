## Environment setup

```r
library(ggplot2)
library(dplyr)
library(scales) # for commas
library(rworldmap) # for map of preprints

setwd('~/code/shithouse/metadata_paper/pub/')

themedarktext = "#707070"
big_fontsize = unit(12, "pt")

country_codes <- data.frame(
  "standard" = c(
    'US',
    'unknown',
    'CN',
    'DK',
    'GB',
    'CA',
    'AU',
    'NL',
    'others'
  ),
  "country" = c(
    'United States',
    'unknown',
    'China',
    'Denmark',
    'United Kingdom',
    'Canada',
    'Australia',
    'Netherlands',
    'others'
  )
)
```

## FIGURE 1

### Panel 1a: World map

Saved as `country_counts.csv`:
```sql
SELECT c.standard, COUNT(DISTINCT t.srs)
FROM roundtwo.tags t
LEFT JOIN roundtwo.standardized_countries c
	ON t.value=c.raw
WHERE t.tag='geo_loc_name'
GROUP BY 1
ORDER BY 2 DESC
```

```r
data <- read.csv('country_counts.csv', header = TRUE, stringsAsFactors = FALSE)
colnames(data) <- c('country','value')
toplot <- joinCountryData2Map(data, joinCode = "ISO2",
  nameJoinColumn = "country",  mapResolution = "coarse", verbose=TRUE)

map <- mapCountryData(toplot, nameColumnToPlot = "value",
        catMethod = "logFixedWidth", numCats=5,
        xlim = NA, ylim = NA, mapRegion = "world",
        colourPalette = "heat", addLegend = FALSE, borderCol = "grey",
        mapTitle = "",
        oceanCol = NA, aspect = 1,
        missingCountryCol = NA, add = FALSE,
        lwd = 0.5)

x <- do.call(
  addMapLegend,
  c(map, legendLabels="all",
    horizontal=FALSE,
    legendShrink = 0.3,
    legendIntervals='page',
    legendMar=5, # move it closer to the continents
    legendWidth=0.7,
    labelFontSize=0.85,
    digits=2
    #tcl=-1.2 # tick mark
  )
)
```


## Panel 1b: Countries over time

Saved as `country_years.csv`:
```sql
WITH top10 AS (
	SELECT c.standard, COUNT(DISTINCT t.srs) AS count
	FROM roundtwo.tags t
	LEFT JOIN roundtwo.standardized_countries c
		ON t.value=c.raw
	WHERE t.tag='geo_loc_name'
	GROUP BY 1
	ORDER BY 2
	LIMIT 11 -- because of unknowns
)

SELECT c.standard, LEFT(s.pubdate,4) AS year,
	COUNT(DISTINCT t.srs) AS count
FROM roundtwo.tags t
LEFT JOIN roundtwo.standardized_countries c
	ON t.value=c.raw
INNER JOIN roundtwo.samples s
	ON t.srs=s.srs
WHERE t.tag='geo_loc_name'
	AND c.standard IN (SELECT standard FROM top10)
GROUP BY 1,2
ORDER BY 1,2
UNION ALL -- then add the "everybody else" category
SELECT 'others' AS standard, LEFT(s.pubdate,4) AS year,
	COUNT(DISTINCT t.srs) AS count
FROM roundtwo.tags t
LEFT JOIN roundtwo.standardized_countries c
	ON t.value=c.raw
INNER JOIN roundtwo.samples s
	ON t.srs=s.srs
WHERE t.tag='geo_loc_name'
	AND c.standard NOT IN (SELECT standard FROM top10)
GROUP BY 1,2
ORDER BY 1,2
```

```r
cyears <- read.csv('country_years.csv', header = TRUE, stringsAsFactors = FALSE)

printable <- cyears %>%
  inner_join(country_codes, by="standard")  %>%
  select(year, country, count)


printable$country <- factor(printable$country,
      levels = c(
        'unknown',
        'others',
        'Netherlands',
        'Australia',
        'Canada',
        'United Kingdom',
        'Denmark',
        'China',
        'United States'
      ))

ggplot(printable, aes(x=year, y=count, fill=country)) +
  geom_bar(stat="identity", color="white") +
  theme_bw() +
  labs(x="Year", y="Samples published") +
  scale_y_continuous(labels=comma) +
  scale_fill_brewer(
    palette='Set1', direction=-1
  ) +
  theme(
    axis.text = element_text(size=big_fontsize, color = themedarktext),
    axis.title = element_text(size=big_fontsize)
  )
```