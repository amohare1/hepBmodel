# Import modules
import random
import csv
import numpy as np

# Variables
total_n = 1000                 # Number of iterations (simulated population)
total_F = 0.5                  # percentage female of total simulated population
max_age = 75*12                  # maximum age (months)
initial_age = 0 * 12
HBV_prev = 1                   # prevalence of HBV in cohort (%)
HBV_F = 0.5                    # % female of HBV population

Sag_loss_prob = 0.004788498    # monthly probability of SAg loss (rate = 1% per year) for those eAg- with low DNA
Sag_loss_DNA = 220             # Maximum threshold of DNA under which SAg can be lost (10^x/100)

eAg_prev = 1                   # prevalence of eAg among those with HBV (%)
eAg_SC_rate = 0.0049046      # probability of eAg sero-clearance (monthly probability)
eAg_min_age = 1 * 12           # minimum age of eAg sero-clearance (months)
eAg_max_age = 55 * 12          # maximum age of eAg sero-clearance (months)
phase2_time_mean = 5*12        # mean length of time in phase 2 (months)
phase2_time_stdev = 6          # standard deviation (months)

p4_min_age = 0 * 12             # minimum age of phase 4 (immune escape)
p3_to_p4_1 = 0.002131059        # probability of transitioning from phase 3 to phase 4 in risk group 1
p3_to_p4_2 = 0.002280729        # probability of transitioning from phase 3 to phase 4 in risk group 2
p3_to_p4_3 = 0.008934845        # probability of transitioning from phase 3 to phase 4 in risk group 3
p3_4_group1 = 30 * 12           # Age distinguishing risk group 1 & 2
p3_4_group2 = 40 * 12           # Age distinguishing risk group 2 & 3

mean_DNA_IT = 800               # average DNA for those in phase 1
stdev_DNA_IT = 35               # standard deviation of DNA for those in phase 1
DNA_p1_amp = 10                 # amplitude of monthly DNA changes (10^x/100), by phase
DNA_p2_amp = 24
DNA_p3_amp = 3                  # DNA_p3_amp should be 0 if DNA stays at a "setpoint' during p3
DNA_p4_amp = 10
DNA_p2_effect = -10
DNA_p4_effect = 0

HCC_prev = 0                    # prevalence of HCC in initial cohort
HCC_risk_months = 12            # Number of prior months that are used to consider DNA level in HCC risk
HCC_risk_DNA = 0                # Average monthly DNA used to determine HCC risk
risk_category_one = 248         # Risk stratification and corresponding HCC probabilities based on DNA level
HCC_cat_one_prob_m = 0.000369925
HCC_cat_one_prob_f = 0.000142295
risk_category_two = 400
HCC_cat_two_prob_m = 0.000469371
HCC_cat_two_prob_f = 0.000180553
risk_category_three = 500
HCC_cat_three_prob_m = 0.000844905
HCC_cat_three_prob_f = 0.000325048
risk_category_four = 600
HCC_cat_four_prob_m = 0.002049459
HCC_cat_five_prob_m = 0.0027246
HCC_cat_four_prob_f = 0.000788751
HCC_cat_five_prob_f = 0.001048803
HCC_cirrhosis_risk_multiplier = 11.7 # HR of HCC for patients with cirrhosis v. no cirrhosis

Cirrhosis_prev = 0              # prevalence of cirrhosis in initial cohort
Cirrhosis_risk_months = 12      # Number of prior months that are used to consider DNA level in cirrhosis risk
Cirrhosis_risk_DNA = 0          # Average monthly DNA used to determine Cirrhosis risk
Cirrhosis_cat_one_prob = 0.000282293 # Risk stratification and corresponding HCC probabilities based on DNA level
Cirrhosis_cat_two_prob = 0.000358186
Cirrhosis_cat_three_prob = 0.000644792
Cirrhosis_cat_four_prob = 0.001564275
Cirrhosis_cat_five_prob = 0.002079751

treatment = False              # treatment incorporated into simulation?

# Translate annual rates into monthly probabilities

# Prep arrays
array_index = max_age
DNA_array = np.zeros((array_index+1, total_n + 1))
event_array = np.zeros((total_n, 8))
phase_DNA_array = np.zeros((max_age * total_n, 3))

# Cohort simulation
total_LM = 0

p1_count = 0            # Number of people remaining in Phase 1 by end of simulation
p3_count = 0            # Number of people in Phase 3 by end of simulation
p4_count = 0            # Number of people in phase 4 by end of simulation
SAg_loss_count = 0      # Number of people losing SAg by end of simulation
eAg_loss_count = 0      # Number of people losing eAg
DNA_lo_count = 0        # Number of people with low (<2000) DNA level
DNA_hi_count = 0        # Number of people with high (>20,000) DNA level


for count in range(0, total_n, 1):

        month = 0
        while month < max_age:  # Individual simulation
                month += 1

                if month == 1:

                        # Initialize age
                        age_months = initial_age

                        # Initialize sex
                        sex = random.randint(1,2)

                        # Assign HBV SAg status and eAg status
                        HBV_sag_rand = random.random()
                        if HBV_sag_rand <= HBV_prev:
                                HBV_sag = True
                                Sag_loss_age = max_age * 2 # initialize SAg_loss_age (age at which SAg is lost)

                                eAg_rand = random.random()
                                if eAg_rand <= eAg_prev:
                                        eAg = True
                                        eAg_age_m = max_age * 2 # initialize age at eAg seroclearance
                                        p4_age = max_age * 2 # initialize age at which phase 4 begins

                                        # Choose phase of illness
                                        phase = 1

                                        # Initialize DNA
                                        DNA = int(np.random.normal(mean_DNA_IT, stdev_DNA_IT, 1))

                                elif eAg_rand > eAg_prev:
                                        eAg = False
                                        eAg_age_m = int(np.random.normal(30, 3.5, 1))

                                        # Choose phase of illness? P3, P4
                                        phase = 3
                                        p4_age = max_age * 2 # initialize age at which phase 4 begins

                                        # Initialize DNA
                                        initial_DNA_rand = random.random()
                                        if initial_DNA_rand <= 0.28:
                                                DNA = int(np.random.normal(220, 10, 1))
                                        elif 0.28 < initial_DNA_rand <= 0.65:
                                                DNA = int(np.random.normal(310, 30, 1))
                                        elif 0.65 < initial_DNA_rand <= 0.85:
                                                DNA = int(np.random.normal(450, 15, 1))
                                        elif 0.85 < initial_DNA_rand <= 1:
                                                DNA = int(np.random.normal(700, 50, 1))

                                # Assign HCC status
                                HCC_rand = random.random()
                                if HCC_rand <= HCC_prev:
                                        HCC = True
                                        # Determine age at which HCC was acquired
                                else:
                                        HCC = False
                                        HCC_age = max_age * 2 # initialize age at which HCC occurs

                                # Assign Cirrhosis status
                                Cirrhosis_rand = random.random()
                                if Cirrhosis_rand <= Cirrhosis_prev:
                                        Cirrhosis = True
                                        # Determine age at which Cirrhosis was acquired
                                else:
                                        Cirrhosis = False
                                        Cirrhosis_age = max_age * 2 # initialize the age at which cirrhosis occurs

                        elif HBV_sag_rand > HBV_prev:
                                HBV_sag = False
                                eAg = False
                                HCC = False
                                Cirrhosis = False

                else:
                        # Age Updater
                        age_months += 1

                        if HBV_sag == True:

                                # eAg updater
                                if phase == 1:
                                        if eAg_min_age < age_months < eAg_max_age:
                                                # Roll for eAg seroclearance
                                                e_SC_rand = random.random()
                                                if e_SC_rand <= eAg_SC_rate:
                                                        # Phase 2 begins
                                                        phase = 2

                                                        # Draw for a phase 2 time
                                                        phase2_time = int(np.random.normal(phase2_time_mean, phase2_time_stdev, 1))

                                                        #Ideally would want to draw phase 2_time from a distribution
                                                        eAg_age_m = age_months + phase2_time

                                                        # DNA updater
                                                        DNA += random.randint(-DNA_p2_amp, DNA_p2_amp) + DNA_p2_effect
                                                        # ALT rule for p2
                                                        # HCC risk for p2
                                                else:
                                                        DNA += random.randint(-DNA_p1_amp, DNA_p1_amp)
                                                        # ALT rule for p1
                                                        # HCC rule for p1
                                        else:
                                                DNA += random.randint(-DNA_p1_amp, DNA_p1_amp)
                                        DNA = max(0, DNA)
                                        # Max DNA? DNA = min(1200, DNA)

                                elif phase == 2:
                                        if age_months < eAg_age_m:
                                                DNA += random.randint(-DNA_p2_amp, DNA_p2_amp) + DNA_p2_effect
                                                # ALT rule for p2
                                                # HCC risk for p2
                                        elif age_months == eAg_age_m:
                                                eAg = False
                                                phase = 3
                                                # Draw for DNA distribution in phase 3
                                                DNA += random.randint(-DNA_p3_amp, DNA_p3_amp)
                                        DNA = max(0, DNA)

                                elif phase == 3:
                                        # SAg updater
                                        if eAg == False: # This if statement is redundant
                                                if DNA < Sag_loss_DNA:
                                                        Sag_loss_rand = random.random()
                                                        if Sag_loss_rand <= Sag_loss_prob:
                                                                HBV_sag = False
                                                                Sag_loss_age = age_months
                                                                phase = 5               # Phase "5" is SAg loss

                                        # DNA stays at its set point at phase 3
                                        if HBV_sag == True:
                                                DNA += random.randint(-DNA_p3_amp, DNA_p3_amp)

                                                # ALT rule for p3: trend towards normal?

                                                # Cirrhosis risk for p3
                                                if sex == 1:
                                                        sum = 0
                                                        Cirrhosis_risk_DNA = 0
                                                        for aa in range(1, Cirrhosis_risk_months+1, 1):
                                                                sum += DNA_array[age_months-aa, count]
                                                        Cirrhosis_risk_DNA = sum/ Cirrhosis_risk_months
                                                        if Cirrhosis_risk_DNA <= risk_category_one:
                                                                Cirrhosis_risk_rand = random.random()
                                                                if Cirrhosis_risk_rand <= Cirrhosis_cat_one_prob:
                                                                        Cirrhosis = True
                                                                        Cirrhosis_age = age_months
                                                        elif risk_category_one < Cirrhosis_risk_DNA <= risk_category_two:
                                                                Cirrhosis_risk_rand = random.random()
                                                                if Cirrhosis_risk_rand <= Cirrhosis_cat_two_prob:
                                                                        Cirrhosis = True
                                                                        Cirrhosis_age = age_months
                                                        elif risk_category_two < Cirrhosis_risk_DNA <= risk_category_three:
                                                                Cirrhosis_risk_rand = random.random()
                                                                if Cirrhosis_risk_rand  <= Cirrhosis_cat_three_prob:
                                                                        Cirrhosis = True
                                                                        Cirrhosis_age = age_months

                                                        # HCC risk for p3
                                                        # Note: re-do 'sum' calculation to allow for the possibility that
                                                        # Cirrhosis_risk_months do not equal HCC_risk_months
                                                        sum = 0
                                                        HCC_risk_DNA = 0
                                                        for aa in range(1, HCC_risk_months+1, 1):
                                                                sum += DNA_array[age_months-aa, count]
                                                        HCC_risk_DNA = sum/ HCC_risk_months
                                                        if HCC_risk_DNA <= risk_category_one:
                                                                HCC_risk_rand = random.random()
                                                                if HCC_risk_rand <= HCC_cat_one_prob_m:
                                                                        HCC = True
                                                                        HCC_age = age_months
                                                        elif risk_category_one < HCC_risk_DNA <= risk_category_two:
                                                                #HCC risk
                                                                HCC_risk_rand = random.random()
                                                                if HCC_risk_rand <= HCC_cat_two_prob_m:
                                                                        HCC = True
                                                                        HCC_age = age_months
                                                        elif risk_category_two < HCC_risk_DNA <= risk_category_three:
                                                                #HCC Risk
                                                                HCC_risk_rand = random.random()
                                                                if HCC_risk_rand <= HCC_cat_three_prob_m:
                                                                        HCC = True
                                                                        HCC_age = age_months
                                                elif sex == 2:
                                                        sum = 0
                                                        Cirrhosis_risk_DNA = 0
                                                        for aa in range(1, Cirrhosis_risk_months + 1, 1):
                                                                sum += DNA_array[age_months - aa, count]
                                                        Cirrhosis_risk_DNA = sum / Cirrhosis_risk_months
                                                        if Cirrhosis_risk_DNA <= risk_category_one:
                                                                Cirrhosis_risk_rand = random.random()
                                                                if Cirrhosis_risk_rand <= Cirrhosis_cat_one_prob:
                                                                        Cirrhosis = True
                                                                        Cirrhosis_age = age_months
                                                        elif risk_category_one < Cirrhosis_risk_DNA <= risk_category_two:
                                                                Cirrhosis_risk_rand = random.random()
                                                                if Cirrhosis_risk_rand <= Cirrhosis_cat_two_prob:
                                                                        Cirrhosis = True
                                                                        Cirrhosis_age = age_months
                                                        elif risk_category_two < Cirrhosis_risk_DNA <= risk_category_three:
                                                                Cirrhosis_risk_rand = random.random()
                                                                if Cirrhosis_risk_rand <= Cirrhosis_cat_three_prob:
                                                                        Cirrhosis = True
                                                                        Cirrhosis_age = age_months

                                                        # HCC risk for p3
                                                        # Note: re-do 'sum' calculation to allow for the possibility that
                                                        # Cirrhosis_risk_months do not equal HCC_risk_months
                                                        sum = 0
                                                        HCC_risk_DNA = 0
                                                        for aa in range(1, HCC_risk_months + 1, 1):
                                                                sum += DNA_array[age_months - aa, count]
                                                        HCC_risk_DNA = sum / HCC_risk_months
                                                        if HCC_risk_DNA <= risk_category_one:
                                                                HCC_risk_rand = random.random()
                                                                if HCC_risk_rand <= HCC_cat_one_prob_f:
                                                                        HCC = True
                                                                        HCC_age = age_months
                                                        elif risk_category_one < HCC_risk_DNA <= risk_category_two:
                                                                # HCC risk
                                                                HCC_risk_rand = random.random()
                                                                if HCC_risk_rand <= HCC_cat_two_prob_f:
                                                                        HCC = True
                                                                        HCC_age = age_months
                                                        elif risk_category_two < HCC_risk_DNA <= risk_category_three:
                                                                # HCC Risk
                                                                HCC_risk_rand = random.random()
                                                                if HCC_risk_rand <= HCC_cat_three_prob_f:
                                                                        HCC = True
                                                                        HCC_age = age_months

                                                # roll dice on p4
                                                # This probability is dependent on eAg_age
                                                if age_months > p4_min_age:
                                                        if eAg_age_m <= p3_4_group1:
                                                                p4_rand = random.random()
                                                                if p4_rand <= p3_to_p4_1:
                                                                        phase = 4
                                                                        p4_age = age_months
                                                        elif p3_4_group1 < eAg_age_m <= p3_4_group2:
                                                                p4_rand = random.random()
                                                                if p4_rand <= p3_to_p4_2:
                                                                        phase = 4
                                                                        p4_age = age_months
                                                        elif p3_4_group2 < eAg_age_m:
                                                                p4_rand = random.random()
                                                                if p4_rand <= p3_to_p4_3:
                                                                        phase = 4
                                                                        p4_age = age_months
                                                DNA = max(0, DNA)

                                else:
                                        DNA += random.randint(-DNA_p4_amp, DNA_p4_amp) + DNA_p4_effect
                                        # ALT rule for p4: probably same rule for p2

                                        # Cirrhosis risk for p3
                                        if sex == 1:

                                                sum = 0
                                                Cirrhosis_risk_DNA = 0
                                                for aa in range(1, Cirrhosis_risk_months + 1, 1):
                                                        sum += DNA_array[age_months - aa, count]
                                                Cirrhosis_risk_DNA = sum / Cirrhosis_risk_months
                                                if Cirrhosis_risk_DNA <= risk_category_four:
                                                        Cirrhosis_risk_rand = random.random()
                                                        if Cirrhosis_risk_rand <= Cirrhosis_cat_four_prob:
                                                                Cirrhosis = True
                                                                Cirrhosis_age = age_months
                                                elif risk_category_four < Cirrhosis_risk_DNA:
                                                        Cirrhosis_risk_rand = random.random()
                                                        if Cirrhosis_risk_rand <= Cirrhosis_cat_five_prob:
                                                                Cirrhosis = True
                                                                Cirrhosis_age = age_months

                                                # HCC risk for p4
                                                sum = 0
                                                HCC_risk_DNA = 0
                                                for aa in range(1, HCC_risk_months + 1, 1):
                                                        sum += DNA_array[age_months - aa, count]

                                                HCC_risk_DNA = sum / HCC_risk_months
                                                if HCC_risk_DNA <= risk_category_four:
                                                        # HCC risk
                                                        HCC_risk_rand = random.random()
                                                        if HCC_risk_rand <= HCC_cat_four_prob_m:
                                                                HCC = True
                                                                HCC_age = age_months
                                                elif risk_category_four < HCC_risk_DNA:
                                                        # HCC risk
                                                        HCC_risk_rand = random.random()
                                                        if HCC_risk_rand <= HCC_cat_five_prob_m:
                                                                HCC = True
                                                                HCC_age = age_months
                                                DNA = max(0, DNA)

                                        elif sex == 2:
                                                sum = 0
                                                Cirrhosis_risk_DNA = 0
                                                for aa in range(1, Cirrhosis_risk_months + 1, 1):
                                                        sum += DNA_array[age_months - aa, count]
                                                Cirrhosis_risk_DNA = sum / Cirrhosis_risk_months
                                                if Cirrhosis_risk_DNA <= risk_category_four:
                                                        Cirrhosis_risk_rand = random.random()
                                                        if Cirrhosis_risk_rand <= Cirrhosis_cat_four_prob:
                                                                Cirrhosis = True
                                                                Cirrhosis_age = age_months
                                                elif risk_category_four < Cirrhosis_risk_DNA:
                                                        Cirrhosis_risk_rand = random.random()
                                                        if Cirrhosis_risk_rand <= Cirrhosis_cat_five_prob:
                                                                Cirrhosis = True
                                                                Cirrhosis_age = age_months

                                                # HCC risk for p4
                                                sum = 0
                                                HCC_risk_DNA = 0
                                                for aa in range(1, HCC_risk_months + 1, 1):
                                                        sum += DNA_array[age_months - aa, count]

                                                HCC_risk_DNA = sum / HCC_risk_months
                                                if HCC_risk_DNA <= risk_category_four:
                                                        # HCC risk
                                                        HCC_risk_rand = random.random()
                                                        if HCC_risk_rand <= HCC_cat_four_prob_f:
                                                                HCC = True
                                                                HCC_age = age_months
                                                elif risk_category_four < HCC_risk_DNA:
                                                        # HCC risk
                                                        HCC_risk_rand = random.random()
                                                        if HCC_risk_rand <= HCC_cat_five_prob_f:
                                                                HCC = True
                                                                HCC_age = age_months
                                                DNA = max(0, DNA)


                        else:
                                #if SAg is false, then compute background mortality ?and non-HBV cirrhosis risk?
                                dummy = True
                DNA = min(1000, DNA)
                # Print output for this month into an array
                DNA_array[age_months, count + 1] = DNA
                phase_DNA_array[count * max_age + age_months, 0] = count
                phase_DNA_array[count * max_age + age_months, 1] = phase
                phase_DNA_array[count * max_age + age_months, 2] = float(DNA/100)

        # Print outputs for this simulated iteration
        #DNA_array[0, count] = count
        event_array[count, 0] = count
        event_array[count, 1] = sex
        event_array[count, 2] = float(eAg_age_m/12)
        event_array[count, 3] = float(Sag_loss_age/12)
        event_array[count, 4] = float(p4_age/12)
        event_array[count, 5] = float(HCC_age/12)
        event_array[count, 6] = float(Cirrhosis_age/12)
        event_array[count, 7] = int(phase)

        # Print a DNA array for persistent eAg positive

        # Print a DNA array for inactive carriers

        # Print a DNA array for eAg negative chronic hepatitis

for row in range(0, max_age, 1):
        DNA_array[row, 0] = row

DNA_array2 = DNA_array.copy()

#for nnn in range(0, count, 1):
#        DNA_array2[0, n] = event_array[n, 5]

# Export to csv
with open('results/result2-DNA.csv', 'w', ) as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for word in DNA_array:
        wr.writerow(word)

with open('results/result1-events.csv', 'w', ) as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    title_row = ['Count', 'Sex', 'eAg Age', 'SAg loss Age', 'P4 Age', 'HCC Age', 'Cirrhosis Age', 'Final Phase']
    wr.writerow(title_row)
    for word in event_array:
        wr.writerow(word)

with open('results/result3-DNA_phase.csv', 'w', ) as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for word in phase_DNA_array:
        wr.writerow(word)