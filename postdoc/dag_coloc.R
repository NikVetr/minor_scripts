library(dagitty)
dag <- dagitty('dag{
    variant [pos="0,1"]
    expr1 [pos="1,0"]
    expr2 [pos="1,1"]
    expr3 [pos="1,2"]
    phenotype [pos="2,1"]
    variant -> expr1 -> phenotype
    variant -> expr2 
    expr3 -> phenotype}'
)

dag <- dagitty('dag{
    scoring_points [pos="-1,0"]
    helping [pos="1,0"]
    money [pos="0,-1"]
    on_team [pos="0,1"]
    on_team <- scoring_points -> money
    on_team <- helping -> money 
    }'
)
plot(dag)
impliedConditionalIndependencies(dag)
adjustmentSets(dag, outcome = "money", exposure = "scoring_points", type = "minimal", effect = "direct")
