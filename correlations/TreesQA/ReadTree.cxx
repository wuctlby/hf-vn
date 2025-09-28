#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <cstring>  // for strcmp
#include <TLeaf.h>
#include <TLeafF.h>
#include <TLeafD.h>
#include <TLeafI.h>
#include <TList.h>
#include <TKey.h>
#include <TObject.h>
#include <TObjArray.h>
#include <map>
#include <string>
#include <TH1D.h>
#include "yaml-cpp/yaml.h"   // Needs yaml-cpp installed

void ReadTree(std::string configFile = "config_read_tree.yaml") {
    // Load YAML config
    YAML::Node config = YAML::LoadFile(configFile.data());

    std::string inputFile = config["InputFile"].as<std::string>();
    TFile *file = TFile::Open(inputFile.data(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    
    TList *folders = file->GetListOfKeys();
    TIter next(folders);
    TObject* folder;
    
    // Remove parentFiles folders
    while ((folder = next())) {
        if (strcmp(folder->GetName(), "parentFiles") == 0) {
            folders->Remove(folder);
        }
    }
    
    // Extract query (assuming 1 key-value pair for simplicity)
    std::map<std::string, std::string> queries;
    for (const auto& treeNode : config["Trees"]) {
        std::string treeName = treeNode["Name"].as<std::string>();
        std::string query = treeNode["Query"].as<std::string>();
        queries[treeName] = query;
    }

    std::string outFile = config["OutputFile"].as<std::string>();
    TFile *outputFile = TFile::Open(outFile.data(), "RECREATE");
    TIter survivedFolders(folders);
    while ((folder = survivedFolders())) {
        for (const auto& [treeName, query] : queries) {

            std::cout << "[" << folder->GetName() << "] Tree: " << treeName << " | Query: " << query << std::endl;

            // Get tree
            TTree *tree = (TTree*)file->Get(Form("%s/%s", folder->GetName(), treeName.c_str()));
            if (!tree) {
                std::cerr << "Tree not found: " << treeName << std::endl;
                return;
            }

            // Apply query: create a new filtered tree
            outputFile->mkdir(Form("%s/%s", folder->GetName(), treeName.c_str()), Form("%s/%s", folder->GetName(), treeName.c_str()), kTRUE);
            outputFile->cd(Form("%s/%s", folder->GetName(), treeName.c_str()));
            TObjArray* leaves = tree->GetListOfLeaves();
            TIter nextLeaf(leaves);
            TLeaf* leaf;
            while ((leaf = (TLeaf*)nextLeaf())) {
                std::string branchName = leaf->GetName();

                double minTreeBranch = tree->GetMinimum(branchName.c_str());
                double maxTreeBranch = tree->GetMaximum(branchName.c_str());
                std::cout << "Processing branch: " << branchName << " with entries in range " << minTreeBranch << " - " << maxTreeBranch << std::endl;

                std::string histName = "h_" + branchName;

                if (leaf->IsA() == TLeafF::Class() || leaf->IsA() == TLeafD::Class() || leaf->IsA() == TLeafI::Class()) {
                    TH1D* h = new TH1D(histName.c_str(), branchName.c_str(), 100, minTreeBranch, maxTreeBranch);
                    if (strcmp(query.c_str(), "") != 0) {
                        tree->Draw((branchName + ">>" + histName).c_str(), "", "goff"); // fill directly from tree
                    } else {
                        tree->Draw((branchName + ">>" + histName).c_str(), query.c_str(), "goff"); // fill directly from tree
                    }
                    h->Write();
                    delete h;
                }
            }
        }
        std::cout << "Finished processing folder: " << folder->GetName() << std::endl;
        std::cout << std::endl;
    }
    outputFile->Close();
    file->Close();
    exit(0);

}
