
// incorporating efficiency stuff.
// abc

functions {
  vector concat_2d_array(array[,] real x) {
    int elements = prod(dims(x));

    vector[elements] y;

    int counter = 1;
    int n_rows = dims(x)[1];
    for (c in 1:dims(x)[2]) {
      y[counter:(counter + n_rows - 1)] = to_vector(x[,c]);
      counter += n_rows;
    }

    return y;
  }


  matrix predict_ordered_logistic(vector alk_par_beta, vector alk_par_eta, matrix mod_mat, int max_a_overall) {
    // print("called: predict_ordered_logistic");
    int num_etas = num_elements(alk_par_eta);

    // matrix[rows(mod_mat),num_etas + 1] pred_mat;
    matrix[rows(mod_mat),max_a_overall] pred_mat;

    // new version
    vector[rows(mod_mat)] mu;

    for (i in 1:rows(mod_mat)) {
      mu[i] = sum(mod_mat[i,1:cols(mod_mat)] .* alk_par_beta');
    }

    // old version
    // vector[rows(mod_mat)] mu = mod_mat * alk_par_beta;



    pred_mat[1:rows(mod_mat),1] = 1 - inv_logit(mu - alk_par_eta[1]);
    for (k in 2:(num_etas)) {
      pred_mat[1:rows(mod_mat),k] = inv_logit(mu - alk_par_eta[k - 1]) - inv_logit(mu - alk_par_eta[k]);
    }
    // old
    pred_mat[1:rows(mod_mat),num_etas+1] = inv_logit(mu - alk_par_eta[num_etas]);

    // new ... would like to go back to the old version because it is vectorized but it seems to be having issues
    // for (i in 1:rows(mod_mat)) {
    //   pred_mat[i,num_etas] = 1 - sum(pred_mat[i,1:(num_etas-1)]);
    // }

    // maybe this will cause problems
    if (max_a_overall > num_etas + 1) {
      for (i in (num_etas + 2):max_a_overall) {
        // vector[rows(mod_mat)] tmp;
        // tmp[1:rows(mod_mat)] = 0.0;
        pred_mat[1:rows(mod_mat),i] = mu * 0;
      }
    }

    return pred_mat;
  }


}

data {
  int<lower=1> N_m;
  int<lower=1> N_l;
  int<lower=1> N_j;
  int<lower=1> N_k;
  int<lower=1> max_t;
  int<lower=1> N_stations;
  int<lower=0> N_knownage_m;
  int<lower=0> N_knownage_l;
  int<lower=1> max_a_overall;
  int<lower=1> N_batches;
  int<lower=1> N_recap_sites;
  int<lower=1> N_last_sites;
  // int<lower=1> N_recap_sites_not_last;
  int<lower=1> N_not_last_sites;
  int<lower=1> N_theta_r;
  int<lower=1> N_theta_c;
  int<lower=1> N_p_r;
  int<lower=1> N_p_c;
  int<lower=1> N_theta_indices_r;
  int<lower=1> N_theta_indices_c;
  int<lower=1> N_theta_par;
  int<lower=1> N_p_indices_r;
  int<lower=1> N_p_indices_c;
  int<lower=1> N_p_par;
  int<lower=1> max_t_p; // need to check what this is
  int N_site_path_length3;
  int N_not_last_sites_rev;
  int N_overall_surv;
  int N_cohort_surv;
  int N_groups;
  int N_overall_surv_r;
  int N_overall_surv_c;
  int N_cohort_surv_r;
  int N_cohort_surv_c;


  array[N_stations] int set_min_a;
  array[N_stations] int set_max_a;
  array[N_j] int max_s_rel;
  array[N_stations] int max_t_recap;
  array[N_j] int next_site;
  array[N_knownage_m] int knownage_m;
  array[N_knownage_l] int knownage_l;
  array[N_recap_sites] int recap_sites;
  array[N_last_sites,1] int last_sites;
  // array[N_recap_sites_not_last,1] int recap_sites_not_last;
  array[N_not_last_sites] int not_last_sites;
  array[N_site_path_length3,1] int site_path_length3;
  array[N_not_last_sites_rev] int not_last_sites_rev;


  array[N_stations] int batches_list_len; // how many entries of batches for this site
  array[N_stations] int site_path_len; // how many entries of sites for this site path

  array[N_m, 8] int m_matrix;
  array[N_l, 6] int l_matrix;
  matrix[N_theta_r,N_theta_c] mod_mat_theta;
  matrix[N_p_r,N_p_c] mod_mat_p;
  array[N_theta_indices_r, N_theta_indices_c] int indices_theta;
  array[N_p_indices_r, N_p_indices_c] int indices_p_obs;

  array[N_overall_surv_r, N_overall_surv_c] int indices_overall_surv;
  array[N_cohort_surv_r, N_cohort_surv_c] int indices_cohort_surv;

  array[N_stations,max(batches_list_len)] int batches_list;
  array[N_stations, max(site_path_len)] int site_path;

  // age stuff

  int<lower=1> mod_mat_a_r;
  int<lower=1> mod_mat_a_c;
  int<lower=1> N_a_parbeta;
  int<lower=1> N_obsageclass;

  array[N_obsageclass] int obsageclass;
  matrix[mod_mat_a_r, mod_mat_a_c-1] mod_mat_a_beta;

  int mod_mat_a_L_r;
  int mod_mat_a_L_c;
  int mod_mat_a_M_r;
  int mod_mat_a_M_c;
  int N_obsageclass_L;
  int N_obsageclass_M;


  matrix[mod_mat_a_L_r, mod_mat_a_L_c-1] mod_mat_a_L;
  matrix[mod_mat_a_M_r, mod_mat_a_M_c-1] mod_mat_a_M;
  array[N_obsageclass_L] int obsageclass_L;
  array[N_obsageclass_M] int obsageclass_M;

  matrix[mod_mat_a_L_r,max_a_overall] fixed_ageclass_l;
  matrix[mod_mat_a_M_r,max_a_overall] fixed_ageclass_m;

  int N_unknownage_l;
  int N_unknownage_m;

  array[N_unknownage_l] int unknownage_l;
  array[N_unknownage_m] int unknownage_m;


  // misc inits
 array[N_groups,N_stations,max_t,
       N_batches,max_a_overall] vector[max_a_overall] inits_theta;



  array[N_groups,N_stations,N_batches] matrix[max_t,max_a_overall] inits_chi_array;

 array[N_groups,N_stations,max_t,
       N_batches,max_a_overall] vector[max_a_overall+1] inits_Theta;



  array[N_groups,N_stations,N_batches,max_t]
      matrix[max_a_overall,max_a_overall] inits_p_obs;


  array[N_groups,N_stations,N_stations,
        max_t,N_batches] matrix[max_t,max_a_overall] inits_lambda_array;

  array[N_stations, max_t] int min_ageclass_mat;

  array[N_stations, max_t] int max_ageclass_mat;


}


parameters {
  // ordered[max_a_overall-1] alk_par_eta; // this may be bad...
  // vector[N_a_parbeta] alk_par_beta;

  row_vector[N_theta_par] theta_params;

  row_vector[N_p_par] p_params;


}



transformed parameters {


}

model {

  for (i in 1:N_theta_par) {
   theta_params ~ normal(0,10); // sigma = 2
  }

  for (i in 1:N_p_par) {
   p_params ~ normal(0,10); // sigma = 2
  }


  // obsageclass ~ ordered_logistic(mod_mat_a_beta * alk_par_beta,alk_par_eta);
  // target += ordered_logistic_lpmf(obsageclass | mod_mat_a_beta * alk_par_beta,alk_par_eta);

   // WAS: theta[a1,a2,s,j,k,b,g] -> theta[g,j,s,b,a1, VECTOR: a2]

 array[N_groups,N_stations,max_t,
       N_batches,max_a_overall] vector[max_a_overall] theta;



   array[N_groups,N_stations,max_t,
       N_batches,max_a_overall] vector[max_a_overall] theta_inv;

    theta = inits_theta;

    array[N_groups,N_stations,
          max_t,N_batches, max_a_overall] vector[max_a_overall + 1] Theta;

  Theta = inits_Theta; // zeros


    // p_obs[g,k,b,t, MATRIX: a1,a2]

  array[N_groups,N_stations,N_batches,max_t]
      matrix[max_a_overall,max_a_overall] p_obs;

  p_obs = inits_p_obs;

  vector[N_theta_indices_r] mu_indices_theta = inv_logit(mod_mat_theta * theta_params');

  for (i in 1:N_theta_indices_r) {

    int a1 = indices_theta[i,1];
    int a2 = indices_theta[i,2];
    int s = indices_theta[i,3];
    int j = indices_theta[i,5];
    // int k = indices_theta[i,6];
    int b = indices_theta[i,7];
    int g = indices_theta[i,8];



    // theta[g,j,s,b,a1, VECTOR: a2]

    theta[g,j,s,b,a1,a2] = mu_indices_theta[i];


  }

  vector[N_p_indices_r] mu_indices_p_obs = inv_logit(mod_mat_p * p_params');

  for (i in 1:N_p_indices_r) {

    int a1 = indices_p_obs[i,1];
    int a2 = indices_p_obs[i,2];
    int t = indices_p_obs[i,4];
    int k = indices_p_obs[i,6];
    int b = indices_p_obs[i,7];
    int g = indices_p_obs[i,8];

    p_obs[g,k,b,t,a1,a2] = mu_indices_p_obs[i];

  }




  for (g in 1:N_groups) {
    for (j in not_last_sites) {
      int k = next_site[j];

      for (s in 1:max_s_rel[j]) {

        for (a1 in 1:set_max_a[j]) {
          for (b in batches_list[j,1:batches_list_len[j]]) {
            Theta[g,j,s,b,a1,1] = theta[g,j,s,b,a1,1];

            for (a2 in 2:max_a_overall) {

              vector[a2-1] tmp_1s = rep_vector(1.0,a2-1);

              vector[a2-1] theta_inv_tmp = tmp_1s - theta[g,j,s,b,a1,1:(a2-1)];


              Theta[g,j,s,b,a1,a2] = theta[g,j,s,b,a1,a2]*prod(theta_inv_tmp);


            }
            vector[max_a_overall] tmp_1s = rep_vector(1.0,max_a_overall);

            vector[max_a_overall] theta_inv_tmp = tmp_1s - theta[g,j,s,b,a1,1:(max_a_overall)];

            Theta[g,j,s,b,a1,max_a_overall + 1] = prod(theta_inv_tmp);


          }

        }
      }
    }
  }




  // lambda_array[g,j,k,t,b, MATRIX: s, a1]
  array[N_groups,N_stations,N_stations,
        max_t,N_batches] matrix[max_t,max_a_overall] lambda_array;

  lambda_array = inits_lambda_array;


  for (g in 1:N_groups) {
    for (j in not_last_sites) {
      int k = next_site[j];
      for (s in 1:max_s_rel[j]) {
        for (b in batches_list[j,1:batches_list_len[j]]) {
          # can start at 1 (bc the theta and thus Theta is already populated with zeros throughout)
          for (a1 in 1:set_max_a[j]) {
            array[2] int tmp_t_opts = {max_t_recap[k],s + (set_max_a[k] - a1)};
            int tmp_max_t = min(tmp_t_opts);
            for (t in s:tmp_max_t) {
              int a2 = a1 + t - s;

              lambda_array[g,j,k,t,b,s,a1] = Theta[g,j,s,b,a1,a2];
            }
          }
        }
      }
    }
  }


  for (g in 1:N_groups) {
    for (j in site_path_length3[1:N_site_path_length3,1]) {
      array[site_path_len[j]] int tmp_ks = site_path[j,1:site_path_len[j]];
      for (n_k in 3:site_path_len[j]) {
        int k = tmp_ks[n_k];
        int k_min1 =tmp_ks[n_k - 1];

        for (b in batches_list[j,1:batches_list_len[j]]) {
          for (s in 1:max_s_rel[j]) {
            for (a1 in 1:set_max_a[j]) {
              array[2] int tmp_vec = {max_t_recap[k],s + (set_max_a[k] - a1)};

              int tmp_max_t = min(tmp_vec);
              for (t in (s):(tmp_max_t)) {
                int a2 = a1 + t - s;
                array[2] int tmp_vec2 = {t - s + a1, set_max_a[k]};
                int tmp_upperage = min(tmp_vec2);

                vector[t - s + 1] tmp_lambda1 = to_vector(lambda_array[g,j,k_min1,s:t,b,s,a1]);


                matrix[1+tmp_upperage-a1,1 + t - s] tmp_p_obs = to_matrix(p_obs[g,k_min1,b,s:t,a1,a1:tmp_upperage]);


                matrix[1+tmp_upperage-a1,1 + t - s] tmp_lambda2 = lambda_array[g,k_min1,k,t,b,s:t,a1:tmp_upperage];


                real tmp_lambda3 = sum(tmp_lambda1 .*
                  ((1 - diagonal(tmp_p_obs)) .* diagonal(tmp_lambda2)));

                lambda_array[g,j,k,t,b,s,a1] = tmp_lambda3;

              }
            }
          }
        }
      }
    }
  }


  array[N_groups,N_stations,N_batches] matrix[max_t,max_a_overall] chi_array;

  chi_array = inits_chi_array;


  for (g in 1:N_groups) {
    for (j in not_last_sites_rev) {
      int k = next_site[j];
      for (s in 1:max_s_rel[j]) {
        for (b in batches_list[j,1:batches_list_len[j]]) {
          # can start at 1 (bc probs are already populated with zeros throughout)
           for (a1 in 1:set_max_a[j]) {
            array[2] int tmp_ts = {max_t_recap[k],s + (set_max_a[k] - a1)};
            int tmp_max_t = min(tmp_ts);

            array[2] int tmp_ages = {set_max_a[k],tmp_max_t - s + a1};
            int tmp_maxage = min(tmp_ages);


            vector[tmp_maxage - a1 + 1] tmp_A = Theta[g,j,s,b,a1,a1:(tmp_maxage)];


            matrix[tmp_maxage - a1 + 1,tmp_max_t - s + 1] tmp_p_obs = to_matrix(p_obs[g,k,b,s:tmp_max_t,a1,a1:tmp_maxage]);


            matrix[tmp_maxage - a1 + 1,tmp_max_t - s + 1] tmp_chi = chi_array[g,k,b,s:tmp_max_t,a1:(tmp_maxage)];


            chi_array[g,j,b,s,a1] = Theta[g,j,s,b,a1,max_a_overall+1] + sum(tmp_A .* (1 - diagonal(tmp_p_obs)) .* diagonal(tmp_chi));


          }
        }
      }
    }
  }





 // matrix[rows(mod_mat_a_L),num_elements(alk_par_eta)+1] pi_L = predict_ordered_logistic(alk_par_beta,
 //                                                                                       alk_par_eta,
 //                                                                                       mod_mat_a_L,
 //                                                                                       max_a_overall);
 //
 // matrix[rows(mod_mat_a_M),num_elements(alk_par_eta)+1] pi_M = predict_ordered_logistic(alk_par_beta,
 //                                                                                       alk_par_eta,
 //                                                                                       mod_mat_a_M,
 //                                                                                       max_a_overall);
 //  // print("pi_M: ",pi_M[1:10,]);

  // print("NLL 1: ",target());

  // print("chi_array: ",chi_array);

  vector[N_l] summand_l;

  for (i in knownage_l) {
    int diff_age = l_matrix[i,2] - l_matrix[i,5];
    int k_age = diff_age + obsageclass_L[i]; //obsageclass_L[i]

     summand_l[i] = log(chi_array[l_matrix[i,4],l_matrix[i,1],l_matrix[i,3],l_matrix[i,2],k_age]) * l_matrix[i,6];

    // summand_l[i] = log(chi_array[l_matrix[i,1],l_matrix[i,2],l_matrix[i,3],l_matrix[i,4],k_age]) * l_matrix[i,6];
  }

  // print("NLL 2: ",sum(summand_l[knownage_l]));

  for (i in unknownage_l) {
    int diff_age = l_matrix[i,2] - l_matrix[i,5];

    // vector[max_a_overall - diff_age] pi_age1 = fixed_ageclass_l[i,(diff_age+1):max_a_overall]';
    vector[max_a_overall - diff_age] pi_age1 = fixed_ageclass_l[i,(1):(max_a_overall - diff_age)]';
    vector[max_a_overall] pi_age2 = append_row(rep_vector(0.0,diff_age),pi_age1)  ./ sum(pi_age1);
    // vector[max_a_overall] pi_age3 = pi_age2 ./ sum(pi_age2);



    // vector[max_a_overall] lik1 = to_vector(chi_array[l_matrix[i,1],l_matrix[i,2],l_matrix[i,3],l_matrix[i,4],1:max_a_overall]);
    // was chi_array[j,s,b,g,a1] -> chi_array[g,j,b, MATRIX: s, a1]
    row_vector[max_a_overall] lik1 = chi_array[l_matrix[i,4],l_matrix[i,1],l_matrix[i,3],l_matrix[i,2],1:max_a_overall];


    // use dot_product
    real lik2 = dot_product(lik1, pi_age2);

    // print("pi: ",pi_age3, "; lik1: ",lik1);

    summand_l[i] = log(lik2) * l_matrix[i,6];

  }

  target += sum(summand_l);
  // print("NLL 3: ",target());

  // print("Chi array: ", chi_array);

  // print("lambda_array: ",lambda_array);


  vector[N_m] summand_m;

  for (i in knownage_m) {
    // diff_age <- (m_matrix[i,"s"] - m_matrix[i,"obs_time"]); diff_age
    // diff_age_t <- (m_matrix[i,"t"] - m_matrix[i,"obs_time"]); diff_age_t
    int diff_age = (m_matrix[i,3] - m_matrix[i,7]);
    int diff_age_t = (m_matrix[i,4] - m_matrix[i,7]);

    // print(max_a-diff_age_t - 1 + 1);
    // print(to_matrix(concat_2d_array(p_obs[1:(max_a-diff_age_t),(1+diff_age_t):max_a,m_matrix[i,4],m_matrix[i,2],m_matrix[i,5]])));

    int k = m_matrix[i,2];
    int s = m_matrix[i,3];
    int t = m_matrix[i,4];

    // vector[max_a_overall-(t - s) - 1 + 1] tmp_p_obs1 = diagonal(to_matrix(
    //     concat_2d_array(p_obs[1:(max_a_overall-(t - s)),(1+ t - s):max_a_overall,m_matrix[i,4],m_matrix[i,2],m_matrix[i,5],m_matrix[i,6]]),
    //     max_a_overall-(t - s) - 1 + 1,max_a_overall-(t - s) - 1 + 1
    //     ));
    // vector[set_max_a[k]-diff_age_t - 1 + 1] tmp_p_obs1 = diagonal(
    //     p_obs[m_matrix[i,6],m_matrix[i,2],m_matrix[i,5],m_matrix[i,4],1:(set_max_a[k]-diff_age_t),(1+diff_age_t):set_max_a[k]]);
    //
    //
    // vector[max_a_overall] tmp_p_obs2 = append_row(tmp_p_obs1,
    //           rep_vector(0.0,max_a_overall - (set_max_a[k]-diff_age_t)));

    int k_age = diff_age + obsageclass_M[i];

    real tmp_p_obs2 = p_obs[m_matrix[i,6],m_matrix[i,2],m_matrix[i,5],m_matrix[i,4],k_age,obsageclass_M[i] + diff_age_t];

    // this may cause some errors. may need to put this in an ifelse statement
    // vector[max_a_overall] tmp_p_obs2 = append_row(rep_vector(0.0,diff_age),
    //           append_row(tmp_p_obs1,rep_vector(0.0,diff_age_t-diff_age)));





    // was lambda_array[j,k,s,t,b,g,a1] -> lambda_array[g,j,k,t,b, MATRIX: s, a1]
    real lik1 =  tmp_p_obs2 * //tmp_p_obs2[k_age + t - s]*
      lambda_array[m_matrix[i,6],m_matrix[i,1],m_matrix[i,2],m_matrix[i,4],m_matrix[i,5],m_matrix[i,3],k_age];
    // real lik1 = tmp_p_obs2[k_age]*
    //   lambda_array[m_matrix[i,1],m_matrix[i,2],m_matrix[i,3],m_matrix[i,4],m_matrix[i,5],m_matrix[i,6],k_age];

    summand_m[i] = log(lik1) * m_matrix[i,8];

    // print("lik1: ",- log(lik1));
    // if (lik1 == 0) {
    //   print("lik == 0");
    //   print(m_matrix[i,]);
    //   print("p = ",tmp_p_obs2);
    //   print("lambda = ",lambda_array[m_matrix[i,6],m_matrix[i,1],m_matrix[i,2],m_matrix[i,4],m_matrix[i,5],m_matrix[i,3],k_age]);
    // }

  }
 // print("NLL 4: ",sum(summand_m[knownage_m]));

 for (i in unknownage_m) {
    int diff_age = (m_matrix[i,3] - m_matrix[i,7]);
    int diff_age_t = (m_matrix[i,4] - m_matrix[i,7]);



    int k = m_matrix[i,2];
    int s = m_matrix[i,3];
    int t = m_matrix[i,4];


    int diff_s_t = t - s;

    vector[max_a_overall-diff_s_t] tmp_p_obs1 = diagonal(
        p_obs[m_matrix[i,6],m_matrix[i,2],m_matrix[i,5],m_matrix[i,4],1:(max_a_overall-diff_s_t),(1+diff_s_t):max_a_overall]);

    vector[max_a_overall] tmp_p_obs2 = append_row(tmp_p_obs1,
              rep_vector(0.0,diff_s_t));

    vector[max_a_overall - diff_age] pi_age1 = fixed_ageclass_m[i,(1):(max_a_overall - diff_age)]';
    vector[max_a_overall] pi_age2 = append_row(rep_vector(0.0,diff_age),pi_age1) ./ sum(pi_age1);
    // vector[max_a_overall] pi_age3 = pi_age2  ./ sum(pi_age2);

    // this may cause some errors. may need to put this in an ifelse statement
    // vector[max_a_overall] tmp_p_obs2 = append_row(rep_vector(0.0,diff_age),
    //           append_row(tmp_p_obs1,rep_vector(0.0,diff_age_t-diff_age)));


    row_vector[max_a_overall] lik1 = tmp_p_obs2' .* lambda_array[m_matrix[i,6],m_matrix[i,1],m_matrix[i,2],m_matrix[i,4],
                                                                      m_matrix[i,5],m_matrix[i,3],1:max_a_overall];


    real lik2 = dot_product(lik1, pi_age2);

    summand_m[i] = log(lik2) * m_matrix[i,8];

    // if (lik2 == 0) {
    //   print("i:", i,"; lik1:", lik1,"; lik2:",lik2,"; tmp_p_obs2:",tmp_p_obs2);
    //   print("pi_age2:",pi_age2);
    //   print("lambda: ",lambda_array[m_matrix[i,1],m_matrix[i,2],m_matrix[i,3],m_matrix[i,4],
    //                                                                   m_matrix[i,5],m_matrix[i,6],1:max_a_overall]);
    // }

 }
 target += sum(summand_m);
 // print("NLL 5: ",target());

}


generated quantities {
  vector[N_overall_surv] overall_surv;
  vector[N_cohort_surv] cohort_surv;

  {
  array[N_groups,N_stations,max_t,
       N_batches,max_a_overall] vector[max_a_overall] theta;

     array[N_groups,N_stations,max_t,
       N_batches,max_a_overall] vector[max_a_overall] theta_inv;




    theta = inits_theta;

    array[N_groups,N_stations,
          max_t,N_batches, max_a_overall] vector[max_a_overall + 1] Theta;

  Theta = inits_Theta; // zeros




  vector[N_theta_indices_r] mu_indices_theta = mod_mat_theta * theta_params';
  for (i in 1:N_theta_indices_r) {

    int a1 = indices_theta[i,1];
    int a2 = indices_theta[i,2];
    int s = indices_theta[i,3];
    int j = indices_theta[i,5];
    // int k = indices_theta[i,6];
    int b = indices_theta[i,7];
    int g = indices_theta[i,8];

    theta[g,j,s,b,a1,a2] += inv_logit(mu_indices_theta[i]);

  }



  for (g in 1:N_groups) {
    for (j in not_last_sites) {
      int k = next_site[j];
      // min_ageclass_mat[j,s]
      for (s in 1:max_s_rel[j]) {
        for (a1 in 1:set_max_a[j]) {
          for (b in batches_list[j,1:batches_list_len[j]]) {
            Theta[g,j,s,b,a1,1] = theta[g,j,s,b,a1,1];
            for (a2 in 2:max_a_overall) {
              // array[a2-1] real tmp_1s;
              // array[a2-1] real tmp_1s;
            vector[a2-1] tmp_1s = rep_vector(1.0,a2-1);


            vector[a2-1] tmp_theta = theta[g,j,s,b,a1,1:(a2-1)];

            vector[a2-1] theta_inv_tmp = tmp_1s - tmp_theta;
            // prod(theta_inv[a1,1:(a2-1),t,j,k,b])


            Theta[g,j,s,b,a1,a2] = theta[g,j,s,b,a1,a2]*prod(theta_inv_tmp);


          }
          vector[max_a_overall] tmp_1s = rep_vector(1.0,max_a_overall);


          vector[max_a_overall] tmp_theta = theta[g,j,s,b,a1,1:(max_a_overall)];

          vector[max_a_overall] theta_inv_tmp = tmp_1s - tmp_theta;
          Theta[g,j,s,b,a1,max_a_overall + 1] = prod(theta_inv_tmp);


        }

      }
    }
   }
  }



  // // not_last_sites <- c(1:n_stations)[-last_sites]
  // int counter = 1;
  // for (g in 1:N_groups) {
  //   for (j in not_last_sites) {
  //     int k = next_site[j];
  //     // min_ageclass_mat[j,s]
  //     for (t in 1:max_s_rel[j]) { // should this be "s"?
  //       for (a1 in 1:set_max_a[j]) {
  //       for (b in batches_list[j,1:batches_list_len[j]]) {
  //
  //
  //         overall_surv[counter] = sum(Theta[a1,a1:set_max_a[k],t,j,k,b,g]);
  //         counter += 1;
  //       }
  //
  //       }
  //     }
  //
  //   }
  // }


  for (i in 1:N_overall_surv_r) {

    int j = indices_overall_surv[i,1];
    int k = indices_overall_surv[i,2];
    int a1 = indices_overall_surv[i,3];
    int s = indices_overall_surv[i,4];
    int b = indices_overall_surv[i,5];
    int g = indices_overall_surv[i,6];

    overall_surv[i] = sum(Theta[g,j,s,b,a1,a1:set_max_a[k]]);

  }

  for (i in 1:N_cohort_surv_r) {

    int a1 = indices_cohort_surv[i,1];
    int a2 = indices_cohort_surv[i,2];
    int s = indices_cohort_surv[i,3];
    // int t = indices_cohort_surv[i,4];
    int j = indices_cohort_surv[i,5];
    int k = indices_cohort_surv[i,6];
    int b = indices_cohort_surv[i,7];
    int g = indices_cohort_surv[i,8];

    cohort_surv[i] = Theta[g,j,s,b,a1,a2];

  }


  } // end thing
}

