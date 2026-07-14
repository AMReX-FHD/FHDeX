- To run the ML model in AMReX need to first install Libtorch (PyTorch C++ library), currently the GNUmakefile is expecting it to be installed in this current directory
-- run this below
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-2.0.1+cpu.zip
mv libtorch libtorch_cpu    

- Run generate_dist_data.py to generate temperature, pressure and theta data from known distributions.

- Run train_ffnn_model.py that uses the generated data above to train a simple feed-forward neural network and save as a TorchScript to be read into the FHDeX code

=== Now we have a model.pt trained and ready to be used by the AMReX code ===

The updated MFsurfchem_functions.cpp has ML model option built in it.

In the inputs_fhd_stag there are flags added at the bottom for stoch_MFsurfchem, use_ml_mfsurfchem and mfsufchem_model_file
