#!/bin/env python
# coding: utf-8

# ## Load required packages

# In[1]:


import numpy as np
import GA_Parallel 
import matplotlib.pyplot
import random
from keras.utils.np_utils import to_categorical
import multiprocessing
from multiprocessing      import Pool, cpu_count
from functools            import partial
import warnings
warnings.filterwarnings('ignore')
import time


# ## Import data

# In[2]:


data = np.load('/home/pinghou/cybertraining2019_team4/research/2019Team4/march_april/march2014_new.npy')


# ## Feature selection by genetic algorithm

# In[ ]:


start = time. time()
sample_size = 10000
data_inputs = data[:sample_size,:-1]
data_outputs = data[:sample_size,-1]
data_outputs = to_categorical(data_outputs, 4)


num_samples = data_inputs.shape[0]
num_feature_elements = data_inputs.shape[1]

train_indices = np.array(random.sample(range(sample_size), int(sample_size*0.7)))
test_indices = np.setdiff1d(range(sample_size),train_indices)
print("Number of training samples: ", train_indices.shape[0])
print("Number of test samples: ", test_indices.shape[0])

"""
Genetic algorithm parameters:
    Population size
    Mating pool size
    Number of mutations
"""
sol_per_pop = 64 # Population size.
num_parents_mating = int(sol_per_pop/2)  # Number of parents inside the mating pool.
num_mutations = 3 # Number of elements to mutate.
num_processor = 32 # Number of processors to parallel
num_generations = 8 
print("Number of population: ",sol_per_pop)
print("Number of processors: ",num_processor)
print("Number of generations: ",num_generations)

pool = multiprocessing.Pool(processes = num_processor)

# Defining the population shape.
pop_shape = (sol_per_pop, num_feature_elements)

# Creating the initial population.
new_population = np.random.randint(low=0, high=2, size=pop_shape)
print(new_population.shape)

best_outputs = []
for generation in range(num_generations):
    print("Generation : ", generation)
    # Measuring the fitness of each chromosome in the population.
    
    if __name__ == '__main__':

        cal_pop_fitness = partial(GA_Parallel.cal_model_fitness, features = data_inputs, 
                                  labels = data_outputs, train_indices = train_indices, 
                                  test_indices = test_indices) 

        fitness = pool.map(cal_pop_fitness, new_population)

    best_outputs.append(np.max(fitness))
    # The best result in the current iteration.
    print("Best result : ", best_outputs[-1])

    # Selecting the best parents in the population for mating.
    parents = GA_Parallel.select_mating_pool(new_population, fitness, num_parents_mating)

    # Generating next generation using crossover.
    offspring_crossover = GA_Parallel.crossover(parents, offspring_size=(pop_shape[0]-parents.shape[0], num_feature_elements))

    # Adding some variations to the offspring using mutation.
    offspring_mutation = GA_Parallel.mutation(offspring_crossover, num_mutations=num_mutations)

    # Creating the new population based on the parents and offspring.
    new_population[0:parents.shape[0], :] = parents
    new_population[parents.shape[0]:, :] = offspring_mutation

# Getting the best solution after iterating finishing all generations.
# At first, the fitness is calculated for each solution in the final generation.
fitness = pool.map(cal_pop_fitness, new_population)
# Then return the index of that solution corresponding to the best fitness.
best_match_idx = np.where(fitness == np.max(fitness))[0]
best_match_idx = best_match_idx[0]

best_solution = new_population[best_match_idx, :]
best_solution_indices = np.where(best_solution == 1)[0]
best_solution_num_elements = best_solution_indices.shape[0]
best_solution_fitness = fitness[best_match_idx]

print("best_match_idx : ", best_match_idx)
print("best_solution : ", best_solution)
print("Selected indices : ", best_solution_indices)
print("Number of selected elements : ", best_solution_num_elements)
print("Best solution fitness : ", best_solution_fitness)
end = time. time()
print(end - start)


# In[4]:


x = range(1,num_generations+1)
matplotlib.pyplot.plot(x,best_outputs)
matplotlib.pyplot.xlabel("Generations")
matplotlib.pyplot.ylabel("Fitness")
matplotlib.pyplot.show()

