# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to implement the AceNeuralEngine class,
which handles the sequence similarity component of the SAT solver.
"""


import numpy as np
import pandas as pd
import torch.nn as nn
import torch
from transformers import BertModel, BertTokenizer
from transformers import AutoTokenizer, AutoModelForMaskedLM


class AceNeuralEngine(nn.Module):
    """
    ACE Neural Engine handles the contextual sequence encoding 
    of the different peptides to generate embeddings. Once the 
    model outputs are calculated the actual represetnation
    can be chosen from the following.


    Representation options:
    - last_hidden_state: last hidden state of the transformer
    - pooler_output: output of the pooler layer
    - cls_embedding: embedding of the [CLS] token
    - mean_pooling: mean pooling of the last hidden state
    - max_pooling: max pooling of the last hidden state

    Representation implementations are based on the HuggingFace Transformers library
    and code from this Kaggle Notebook: https://www.kaggle.com/code/rhtsingh/utilizing-transformer-representations-efficiently
    """

    def __init__(self, base_model, tokenizer, device=None):
        super(AceNeuralEngine, self).__init__()
        self.model = base_model
        self.tokenizer = tokenizer
        self.device = device if device is not None else torch.device('cpu')
        self.model.to(self.device)

    def forward(self, inputs, representation='last_hidden_state'):
        """
        Implements the forward function of the nerural engine to get the output of the model
        before grabbing the embedding as the speccified representation.
        """
        
        # 1. Move the inputs to the correct device
        ### inputs: [batch_size, max_seq_len]
        if isinstance(inputs, list):
            inputs = self.tokenizer(inputs, padding=True, return_tensors='pt')
        inputs = inputs.to(self.device)
        attention_mask = inputs['attention_mask']
        model_outputs = self.model(**inputs)

        # 2. Get the correct transformer representation
        if representation == 'last_hidden_state':    
            representation = model_outputs.hidden_states[-1]
        
        elif representation == 'pooler_output':
            if hasattr(model_outputs, 'pooler_output'):
                # repesentation: [batch_size, pooler_dim]
                representation = model_outputs.pooler_output
            else:
                raise ValueError("Base model does not have a pooler_output")
        
        elif representation == 'cls_embedding':
            # repesentation: [batch_size, hidden_size]
            representation = model_outputs.hidden_states[-1][:, 0, :]
    
        elif representation=='mean_pooling':
            # repesentation: [batch_size, hidden_size]
            input_mask_expanded = attention_mask.unsqueeze(-1).expand(model_outputs.hidden_states[-1].size()).float()
            sum_embeddings = torch.sum(model_outputs.hidden_states[-1] * input_mask_expanded, 1)
            sum_mask = input_mask_expanded.sum(1)
            sum_mask = torch.clamp(sum_mask, min=1e-9)
            mean_embeddings = sum_embeddings / sum_mask
            representation = model_outputs.hidden_states[-1].mean(dim=1)

        elif representation=='max_pooling':
            representation = torch.max(model_outputs.hidden_states[-1], 1)[0]
                
        elif representation == 'mean_max_pooling':
            # repesentation: [batch_size, 2*hidden_size]
            input_mask_expanded = attention_mask.unsqueeze(-1).expand(model_outputs.hidden_states[-1].size()).float()
            sum_embeddings = torch.sum(model_outputs.hidden_states[-1] * input_mask_expanded, 1)
            sum_mask = input_mask_expanded.sum(1)
            sum_mask = torch.clamp(sum_mask, min=1e-9)
            mean_embeddings = sum_embeddings / sum_mask
            max_embeddings = torch.max(model_outputs.hidden_states[-1], 1)[0]
            representation = torch.cat((mean_embeddings, max_embeddings), 1)
        
        elif representation == 'concatenate_pooling':
            # Strategy: concatenate the last four layers and then pool
            # Alternative strategies: concatenate the last two layers and then pool, 
            # intercalate different layers and pool, etc.

            # repesentation: [batch_size, num_final_layers*hidden_size] 
            queried_hidden_states = torch.stack(model_outputs.hidden_states[-4:])
            representation = torch.cat(tuple(queried_hidden_states), dim=-1)

        #TODO: implement other pooling strategies (e.g. weighted_layer_pooling, attention_pooling)
        # elif self.representation == 'weighted_layer_pooling':
        # elif self.representation == 'attention_pooling':
            
        return representation

    def load_weights(self, weights_path):
        """Load weights from a file"""
        self.load_state_dict(torch.load(weights_path, map_location=self.device))

    def save_weights(self, weights_path):
        """Save weights to a file"""
        torch.save(self.model.state_dict(), weights_path)

    def freeze(self):
        """Freeze the model"""
        for param in self.model.parameters():
            param.requires_grad = False

    def unfreeze(self):
        """Unfreeze the model"""
        for param in self.model.parameters():
            param.requires_grad = True

    @staticmethod
    def cosine_similarity(emb1, emb2):
        """Cosine similarity between two vectors"""
        a = emb1.reshape(-1)
        b = emb2.reshape(-1)
        return np.dot(a,b.T)/(np.linalg.norm(a)*np.linalg.norm(b))

    @staticmethod
    def euclidean_similarity(emb1, emb2):
        """Euclidean similarity between two vectors"""
        a = emb1.reshape(-1)
        b = emb2.reshape(-1)
        return 1 - np.linalg.norm(a-b)/(np.linalg.norm(a)+np.linalg.norm(b))

    def embed_sequences(self, sequences, representation='last_hidden_state'):
        """Calculate embeddings for a list of sequences"""
        # Tokenize sequences
        tokenized = self.tokenizer(list(sequences), padding=True, return_tensors='pt')

        # Get embeddings
        with torch.no_grad():
            output = self.forward(tokenized, representation=representation)
            embeddings = output.cpu().numpy()

        assert len(embeddings) == len(sequences)
        return embeddings

    def find_paired_peptides(self, peptide_ids, peptide_sequences, representation='last_hidden_state', sim_fxn='euclidean', threshold=0.65):
        embeddings = self.embed_sequences(peptide_sequences, representation=representation)
        paired_peptide_ids = []
        for i in range(len(peptide_ids)):
            for j in range(len(peptide_ids)):
                if i == j:
                    continue
                if sim_fxn == 'euclidean':
                    metric =self.euclidean_similarity(embeddings[i], embeddings[j])
                    if metric >= threshold:
                        if (peptide_ids[j], peptide_ids[i], metric) in paired_peptide_ids:
                            continue
                        else:
                            paired_peptide_ids.append((peptide_ids[i], peptide_ids[j], metric))
                elif sim_fxn == 'cosine':
                    metric = self.cosine_similarity(embeddings[i], embeddings[j])
                    if metric >= threshold:
                        if (peptide_ids[j], peptide_ids[i], metric) in paired_peptide_ids:
                            continue
                        else:
                            paired_peptide_ids.append((peptide_ids[i], peptide_ids[j], metric))
                else:
                    raise ValueError("Similarity function must be 'euclidean' or 'cosine'")
        return list(paired_peptide_ids)

    @staticmethod
    def to_paired_peptide_df(paired_peptide_triples):
        """Convert a list of paired peptides to a dataframe"""
        return pd.DataFrame(paired_peptide_triples, columns=['peptide_id_1', 'peptide_id_2', 'similarity'])
