
compute_joint_distribution = function(edit_table, site_number_a, site_number_b){

  assertthat::is.number(site_number_a)
  assertthat::is.number(site_number_b)

  site_a = paste0("Site", site_number_a)
  site_b = paste0("Site", site_number_b)

  joint_distribution = as.data.frame(table(edit_table[, c(site_a, site_b)]))

  return(joint_distribution)
}

compute_marginal_distribution = function(edit_table, site_number){
  assertthat::is.number(site_number)

  site = paste0("Site", site_number)
  marginal = as.data.frame(table(edit_table[, site]))

  return(marginal)
}

conditional_a_given_b = function(edit_table, site_number_a, site_number_b){

  assertthat::is.number(site_number_a)
  assertthat::is.number(site_number_b)

  joint = compute_joint_distribution(edit_table, site_number_b, site_number_a)
  marginal = compute_marginal_distribution(edit_table, site_number_b)

  conditional = joint[, 1:2]
  conditional = cbind(conditional, joint$Freq / marginal$Freq)
  colnames(conditional)[3] = "Conditional"

  return(conditional)

}
# small test
#site1 = "ACCGGA"
#site2 = "TCGGGA"
#p_s2_given_s1 = p_1_and_2[which(p_1_and_2$Site1 == site1 & p_1_and_2$Site2 == site2), "Freq"] / p_1[p_1$Var == site1, "Freq"]
#assertthat::are_equal(p_s2_given_s1, p_2_given_1[which(p_2_given_1$Site1 == site1 & p_2_given_1$Site2 == site2), "Prob"])


plot_conditional_a_given_b = function(conditional, site_a, site_b){

  g_cond = ggplot(data = conditional,
                     aes_string(x=site_a, y="Conditional"))+
    facet_grid(site_b ) +
    geom_col()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90))

  return(g_cond)
}
