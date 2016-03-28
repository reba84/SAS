*Program used to estimate p-values for Bayesian statistics;
options ls=80;
proc iml;

post_data=0.6014; * observed posterior for PGM-2 progeny ratios from fresh gels;
sample_size=127;  *sample size for PGM-2 fresh gels;

*expected distribution from tetrasomic hypothesis;
tetra_a=(1/12);
tetra_b=(5/12);
tetra_c=(5/12);
tetra_d=(1/12);

*expected distribution from disomic 1 hypothesis;
di_a=(1/8);
di_b=(3/8);
di_c=(3/8);
di_d=(1/8);

*expected distribution from disomic 2 hypothesis;
di2_a=(0);
di2_b=(1/2);
di2_c=(1/2);
di2_d=(0);

*Set initial values for loops;
pr=0;
pr_wt=0;
sample=0;
step=0;
seed=0;
iteration=10000;    
post_count=0;
post_count_wt=0;

*start while loop for 10000 iterations to calculate probability of a posterior greater than posterior_data;
    do while (step<iteration);
    step=step+1;
*setting loop values for calculations inside 10000 loop;
    pr_tetra=0;
    pr_di=0;
    post_di=0;
    n=0;
    seed=seed+1;
    a_count=0;
    b_count=0;
    c_count=0;
    d_count=0;


    *start while loop to generate random data for disomic distribution;
        do while (n<sample_size);
            n=n+1;
  *random number generator;
            sample = ranuni(seed);

            sample = ceil(sample*8);

    *scoring of genotypes for random data;
            if (sample=1) then a_count=(a_count+1);
            if (2<=sample & sample<=4) then b_count=(b_count+1);
            if (5<=sample & sample<=7) then c_count=(c_count+1);
            if (sample=8) then d_count=(d_count+1);  
        end; *ending the do while loop;

  *calculate p(data|hypothesis tetrasomic);
      fact_n=fact(n)/(fact(a_count)*fact(b_count)*fact(c_count)*fact(d_count));
        pr_tetra=fact_n*(
                          (tetra_a**(a_count))*
                          (tetra_b**(b_count))*
                          (tetra_c**(c_count))*
                          (tetra_d**(d_count))
                         );
*calculate p(data|hypothesis disomic 1);
        pr_di=fact_n*(
                       (di_a**(a_count))*
                       (di_b**(b_count))*
                       (di_c**(c_count))*
                       (di_d**(d_count)) 
                      );
* calculate p(data|hypothesis disomic 2);
        pr_di2=fact_n*(
                       (di2_a**(a_count))*
                       (di2_b**(b_count))*
                       (di2_c**(c_count))*
                       (di2_d**(d_count))
                     );
    
*sum of pr_di and pr_tetra;
        sum_pr=pr_tetra+pr_di+pr_di2;
                                sum_wt_pr=(3/6)*pr_tetra+(2/6)*pr_di+(1/6)*pr_di2;
*set numerator values for likelihood formula;
        numerator=pr_di;
                                numerator_wt=(2/6)*pr_di;
*calculate disomic posterior;
    post_di=numerator/sum_pr;
                                post_di_wt=numerator_wt/sum_wt_pr;

*count number of posteriors greater than the one I obtained from my experimental 
observations(=post_data);

        if post_di < post_data then do;
            post_count = post_count+1;
        end;  *end post_count equal loop;

                                if post_di_wt < post_data then do;
                                    post_count_wt = post_count_wt+1;
                                end; *end post_count weighted loop;

    end; *end while statement for 100000 loop;    

*calculate probability of getting a posterior greater than post_data from the random sample in 10000 iterations;
pr=post_count/iteration;
pr_wt=post_count_wt/iteration;
print pr, pr_wt;

            quit;
                    run;
