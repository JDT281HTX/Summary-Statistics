import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st 

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path) 

# Combine the data into a single DataFrame
combined_data = pd.merge(mouse_metadata, study_results, on="Mouse ID")

# Display the data table for preview 
print(combined_data) 

# Checking the number of mice.
total_mice = combined_data.shape[0] 
print(f"Total number of mice: {total_mice}") 

# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
# Optional: Get all the data for the duplicate mouse ID.
duplicate_mice = combined_data[combined_data.duplicated(subset=['Mouse ID', 'Timepoint'], keep=False)]  
print(duplicate_mice) 

# Create a clean DataFrame by dropping the duplicate mouse by its ID.
duplicate_mouse_ids = duplicate_mice['Mouse ID'].unique()  
clean_data = combined_data[~combined_data['Mouse ID'].isin(duplicate_mouse_ids)] 
print(clean_data) 

# Checking the number of mice in the clean DataFrame.
unique_mice = clean_data.shape[0] 
print(f"Unique number of mice: {unique_mice}") 

# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Use groupby and summary statistical methods to calculate the following properties of each drug regimen:
# mean, median, variance, standard deviation, and SEM of the tumor volume.
# Assemble the resulting series into a single summary DataFrame. 
regimen_stats = clean_data.groupby('Drug Regimen').agg(
    mean_tumor_volume=('Tumor Volume (mm3)', 'mean'),
    median_tumor_volume=('Tumor Volume (mm3)', 'median'),
    variance_tumor_volume=('Tumor Volume (mm3)', 'var'),
    std_tumor_volume=('Tumor Volume (mm3)', 'std'),
    sem_tumor_volume=('Tumor Volume (mm3)', 'sem')
) 
print(regimen_stats) 

# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas. 
regimen_counts = clean_data['Drug Regimen'].value_counts() 
regimen_counts.plot(kind='bar', figsize=(10,6), color='skyblue') 
plt.title('Total Number of Mouse ID/Timepoints for Each Drug Regimen')
plt.xlabel('Drug Regimen')
plt.ylabel('Number of Mouse ID/Timepoints') 
plt.show() 

# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot. 
plt.figure(figsize=(10,6))
plt.bar(regimen_counts.index, regimen_counts.values, color='skyblue') 
plt.title('Total Number of Mouse ID/Timepoints for Each Drug Regimen', fontsize=14)
plt.xlabel('Drug Regimen', fontsize=12)
plt.ylabel('Number of Mouse ID/Timepoints', fontsize=12)
plt.xticks(rotation=45) 
plt.tight_layout()
plt.show() 

# Generate a pie chart, using Pandas, showing the distribution of unique female versus male mice used in the study 
# Get the unique mice with their gender
# Make the pie chart
gender_mice = clean_data.drop_duplicates(subset=['Mouse ID']) 
gender_counts = gender_mice['Sex'].value_counts() 
gender_counts.plot(kind='pie', autopct='%1.1f%%', figsize=(6,6), startangle=90, colors=['lightcoral', 'skyblue'])
plt.title('Distribution of Unique Female vs Male Mice') 
plt.ylabel('')  
plt.show()

# Generate a pie chart, using pyplot, showing the distribution of unique female versus male mice used in the study
# Get the unique mice with their gender
# Make the pie chart 
plt.figure(figsize=(6,6)) 
plt.pie(gender_counts.values, labels=gender_counts.index, autopct='%1.1f%%', startangle=90, colors=['lightcoral', 'skyblue']) 
plt.title('Distribution of Unique Female vs Male Mice') 
plt.show() 

# Calculate the final tumor volume of each mouse across four of the treatment regimens:
# Capomulin, Ramicane, Infubinol, and Ceftamin
# Start by getting the last (greatest) timepoint for each mouse
# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint 
last_timepoints = clean_data.groupby('Mouse ID').agg({'Timepoint': 'max'}).reset_index()
merged_last_timepoints = pd.merge(last_timepoints, clean_data, on=['Mouse ID', 'Timepoint']) 
final_tumor_volumes = merged_last_timepoints[merged_last_timepoints['Drug Regimen'].isin(['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin'])] 
print(final_tumor_volumes) 

# Put treatments into a list for for loop (and later for plot labels)
treatments = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_data = []

# Calculate the IQR and quantitatively determine if there are any potential outliers.
for treatment in treatments:

    # Locate the rows which contain mice on each drug and get the tumor volumes
    tumor_volumes = final_tumor_volumes[final_tumor_volumes['Drug Regimen'] == treatment]['Tumor Volume (mm3)']

    # add subset
    tumor_vol_data.append(tumor_volumes)

    # Determine outliers using upper and lower bounds 
    Q1 = tumor_volumes.quantile(0.25)
    Q3 = tumor_volumes.quantile(0.75)
    IQR = Q3 - Q1 
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR  

# Generate a box plot that shows the distribution of the tumor volume for each treatment group. 
plt.figure(figsize=(10, 6))
plt.boxplot(tumor_vol_data, labels=treatments) 
plt.title('Distribution of Tumor Volume by Treatment Group', fontsize=14)
plt.xlabel('Drug Regimen', fontsize=12)
plt.ylabel('Tumor Volume (mm3)', fontsize=12) 
plt.grid(axis='y')
plt.show() 

# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin 
mouse_id = 'x401' 
capomulin_data = clean_data[(clean_data['Mouse ID'] == mouse_id) & (clean_data['Drug Regimen'] == 'Capomulin')] 
plt.figure(figsize=(10, 6))
plt.plot(capomulin_data['Timepoint'], capomulin_data['Tumor Volume (mm3)'], marker='o', color='blue', label='Tumor Volume') 
plt.title(f'Tumor Volume vs. Time Point for Mouse {mouse_id} Treated with Capomulin', fontsize=14)
plt.xlabel('Time Point (days)', fontsize=12)
plt.ylabel('Tumor Volume (mm3)', fontsize=12)
plt.xticks(capomulin_data['Timepoint'])  
plt.grid()
plt.legend()
plt.show() 

# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen 
capomulin_data2 = clean_data[clean_data['Drug Regimen'] == 'Capomulin'] 
avg_tumor_volume = capomulin_data2.groupby('Weight (g)').agg({'Tumor Volume (mm3)': 'mean'}).reset_index() 
plt.figure(figsize=(10, 6))
plt.scatter(avg_tumor_volume['Weight (g)'], avg_tumor_volume['Tumor Volume (mm3)'], color='orange', marker='o') 
plt.title('Mouse Weight vs. Average Tumor Volume for Capomulin Regimen', fontsize=14)
plt.xlabel('Weight (g)', fontsize=12)
plt.ylabel('Average Tumor Volume (mm3)', fontsize=12)
plt.grid() 
plt.show() 

# Calculate the correlation coefficient and a linear regression model
# for mouse weight and average observed tumor volume for the entire Capomulin regimen 
capomulin_data3 = clean_data[clean_data['Drug Regimen'] == 'Capomulin'] 
avg_tumor_volume = capomulin_data3.groupby('Weight (g)').agg({'Tumor Volume (mm3)': 'mean'}).reset_index() 
correlation_coefficient, _ = st.pearsonr(avg_tumor_volume['Weight (g)'], avg_tumor_volume['Tumor Volume (mm3)'])
slope, intercept, r_value, p_value, std_err = st.linregress(avg_tumor_volume['Weight (g)'], avg_tumor_volume['Tumor Volume (mm3)'])
x_values = avg_tumor_volume['Weight (g)']
regression_line = slope * x_values + intercept

plt.figure(figsize=(10, 6))
plt.scatter(avg_tumor_volume['Weight (g)'], avg_tumor_volume['Tumor Volume (mm3)'], color='orange', marker='o', label='Average Tumor Volume')
plt.plot(x_values, regression_line, color='blue', label='Regression Line') 
plt.title('Mouse Weight vs. Average Tumor Volume for Capomulin Regimen', fontsize=14)
plt.xlabel('Weight (g)', fontsize=12)
plt.ylabel('Average Tumor Volume (mm3)', fontsize=12)
plt.legend()
plt.grid() 
plt.show() 
print(f'Correlation coefficient: {correlation_coefficient}')