import _ from 'lodash'
import chroma from 'chroma-js'
// import $ from 'jquery'
// import 'select2'
import MolArt from 'molart'

window.chroma = chroma

import '@fortawesome/fontawesome-free/css/all.css'
import 'simple-line-icons/css/simple-line-icons.css'
// import 'select2/dist/css/select2.min.css'
import '../css/ndd.css'

import esNddData from '../../data/v2/list.json'

const METHOD = {
    PDB: {
        id: 'PDB',
        label: 'Monomeric structure from Protein Data Bank'
    },
    PDB_MULTI: {
        id: 'PDB_MULTI',
        label: 'Complex structure from Protein Data Bank'
    },    
    ALPHAFOLD: {
        id: "ALPHAFOLD",
        label: 'Predicted structures from AlphaFold'
    }
};

const dataDir = 'data/'

const containerId = 'molartContainer';

const annotationFields = {
    HotSpot3D: {
        name: "Essential3D",
        label: 'Essential sites',
        tooltip: 'Essential sites',
        category: "Structural based annotations"
        , categorical: true
        , derived: false
        ,color: '#FF5135'
        , variant : false
    },
    paraz3dscore: {
        name: "pConservation3D",
        label: 'Paralog conserved sites',
        tooltip: 'Paralog conserved sites',
        category: "Structural based annotations"
        , categorical: false
        , derived: false
        ,color: '#ED7D31'
        , variant : false
    },
    mtr3dscore: {
        name: "mIntolerance3D",
        label: 'Missense constraint sites',
        tooltip: 'Missense constraint sites',
        category: "Structural based annotations"
        , categorical: false
        , derived: false
        ,color: '#9C45E1'
        , variant : false
    },
    PER_3D: {
        name: "PVE3D",
        label: 'pvEnriched3D',
        tooltip: 'pvEnriched3D',
        category: "Structural based annotations"
        ,color: '#E858EC'
        , categorical: true
        , derived: false
        , variant : false
    },
    denovoDB: {
        name: "denovoDB",
        label: 'denovoDB',
        tooltip: 'denovoDB',
        category: "Variants"
        , categorical: true
        , derived: false
        ,color: '#945200'
        , variant : true
    },
    ASD: {
        name: "ASD",
        label: 'ASD (Satterstrom et al)',
        tooltip: 'ASD (Satterstrom et al)',
        category: "Variants"
        , categorical: true
        , derived: false
        ,color: '#16FCFE'
        , variant : true
    },
    DD: {
        name: "DD",
        label: 'DD (Kaplanis et al)',
        tooltip: 'DD (Kaplanis et al)',
        category: "Variants"
        , categorical: true
        , derived: false
        ,color: '#FFFB00'
        , variant : true
    },    
    DEE: {
        name: "EPI25_DEE",
        label: 'DEE (Epi25 Collaborative)',
        tooltip: 'DEE (Epi25 Collaborative)',
        category: "Variants"
        , categorical: true
        , derived: false
        ,color: '#FF8AD8'
        , variant : true
    },
    UKBiobank: {
        name: "UKBiobank",
        label: 'UKBiobank',
        tooltip: 'UKBiobank',        
        category: "Variants"
        , categorical: true
        , derived: false
        ,color: '#5D9B48'
        , variant : true
    },
    pscount: {
        name: "ClinVar_HGMD_count",
        label: 'ClinVar/HGMD (pathogenic)',
        tooltip: 'Pathogenic variants',
        category: "Variants"
        , categorical: false
        , derived: false
        ,color: '#EC0175'
        , variant : true
    },
    gscount: {
        name: "gnomAD_count",
        label: 'gnomAD (population)',
        tooltip: 'Neutral variants',
        category: "Variants"
        , categorical: false
        , derived: false
        ,color: '#0000DE'
        , variant : true
    },
}

const customConfig = generateCustomConfig();

class DataUris {
    constructor() {
        this.annotations = {}; //one for every chain in the structure file, there is one annotations file
        this.structure = undefined;
    }
}

class SSMapping {
    constructor() {
        this.mapping = [];
    }

    get length() {
        return this.mapping.length;
    }

    getRegion(ix) {
        return this.mapping[ix];
    }

    addRegion(seqBegin, seqEnd, strBegin, strEnd) {
        this.mapping.push({
            seq: {
                begin: seqBegin,
                end: seqEnd
            },
            str: {
                begin: strBegin,
                end: strEnd
            }
        })
    }

    minIx(type){
        let minIx = Number.MAX_SAFE_INTEGER;
        this.mapping.forEach(d => {
            if (d[type].begin < minIx){
                minIx = d[type].begin;
            }
        });

        return minIx;
    }

    maxIx(type){
        let maxIx = Number.MIN_SAFE_INTEGER;
        this.mapping.forEach(d => {
            if (d[type].end > maxIx){
                maxIx = d[type].end;
            }
        });

        return maxIx;
    }

    get minIxSeq(){
        return this.minIx("seq");
    }

    get minIxStr(){
        return this.minIx("str");
    }

    get maxIxSeq(){
        return this.maxIx("seq");
    }

    get maxIxStr(){
        return this.maxIx("str");
    }
}

let geneSelection;
let methodSelection;
let sequenceSelection;
let structureSelection;

const assert = function(condition, message) {
    if (!condition)
        throw Error('Assert failed: ' + (message || ''));
};


$( document ).ready(function() {

    geneSelection = $('#geneSelection');
    methodSelection = $('#methodSelection');
    sequenceSelection = $('#sequenceSelection');
    structureSelection = $('#structureSelection');

    geneSelection.select2({
        placeholder: "Select a gene",
        allowClear: false
    });
    methodSelection.select2({
        placeholder: "Select a method",
        allowClear: false
    });
    // sequenceSelection.select2({
    //     placeholder: "Select a sequence",
    //     allowClear: false
    // });
    structureSelection.select2({
        placeholder: "Select a structure",
        allowClear: false
    });
    // return loadData().then(data => {
    //     console.log('data', data);
    //     populateGeneList(data);
    //     $('#btnShow').on("click", () => btnShowOnClick(data))
    // });

    $("#cntNDD").text(Object.keys(esNddData).length);
    populateGeneList(esNddData);
    $('#btnShow').on("click", () => btnShowOnClick(esNddData))
});

function removeRows(tab, ixs) {

    ixs.sort((a,b) => parseInt(b)-parseInt(a)); //no idea why it does not work without parseInt as the ixs should indeed be numbers
    
    ixs.forEach(ix => {
        Object.values(tab).forEach(col => col.splice(ix, 1))
    });
}

function parseTSV(data, filters) {

    const tab = {};
    let colNames = [];
    let colNameIx = {};
    data.split('\n').forEach((line, ix) => {
        const lineTrim = line.trim();
        if (lineTrim === "") return;
        const sLine = lineTrim.split('\t');
        if (sLine.length === 0) return;
        if (ix === 0) {
            colNames = sLine;
            colNames.forEach((d, ix) => {
                tab[d] = [];
                colNameIx[d] = ix;
            });
        } else {
            let matches = true;
            if (filters !== undefined) {
                Object.keys(filters).forEach(k => {
                    if (filters[k] !== sLine[colNameIx[k]]) {
                        matches = false;
                    }
                })
            }
            if (matches) sLine.forEach((d, ix) => tab[colNames[ix]].push(d));
        }
    });

    //in v2 there are duplicities caused by merging procedure upstream, so these need to be removed otherwise it messes with the following procedure
    let seqIxs = tab['Uniprot_position'].map(d => parseInt(d));
    assert(seqIxs.length > 0)
    const dupIxs = [];
    for (let i = 1; i < seqIxs.length; i++) {
        if (seqIxs[i] == seqIxs[i-1]) dupIxs.push(i)
    }
    removeRows(tab, dupIxs);

    return tab;
}

// Get sequence-structure mapping
function getSSMapping(params) {
    return $.get(params.mappingUri).then(data => {
        const tab = parseTSV(data, params.filters);
        let seqIxs = tab['Uniprot_position'].map(d => parseInt(d));
        let strIxs = tab['Position_in_structure'].map(d => parseInt(d));
        
        assert(seqIxs.length === strIxs.length, "Sequence and structure lengths should match");

        let ixLastBreak = 0;
        let ixCurrent = 1;

        const mapping = new SSMapping();

        while (ixCurrent < seqIxs.length) {
            if (seqIxs[ixCurrent - 1] !== seqIxs[ixCurrent] - 1 || strIxs[ixCurrent - 1] !== strIxs[ixCurrent] - 1) {
                mapping.addRegion(seqIxs[ixLastBreak], seqIxs[ixCurrent - 1], strIxs[ixLastBreak], strIxs[ixCurrent - 1]);
                ixLastBreak = ixCurrent;
            }
            ixCurrent += 1;
        }
        mapping.addRegion(seqIxs[ixLastBreak], seqIxs[ixCurrent - 1], strIxs[ixLastBreak], strIxs[ixCurrent - 1]);

        const seqStart = mapping.minIxSeq;
        const seqEnd = mapping.maxIxSeq;
        const strStart = mapping.minIxStr;
        const strEnd = mapping.maxIxStr;

        const mappingJson = {
            id: params.pdbId,
            chainId: params.chainId,
            structure: {
                format: 'PDB', //valid parameters are PDB and mmCIF
                uri: params.structureUri
            },
            start: 1, // where the structure begins with respect to the full molecule sequence
            end: Math.max(seqEnd, strEnd), // where the structure ends with respect to the full molecule sequence
            seqStart: seqStart, // where the structure begins with respect to the sequence (the sequence does not have to covert the full molecule, therefore seqStart does not have to be 1)
            seqEnd: seqEnd, // where the structure ends with respect to the sequence
            coverage: []
        };

        for (let i = 0; i < mapping.length; ++i) {
            const m = mapping.getRegion(i);
            mappingJson.coverage.push({
                start: {
                    residue_number: m.seq.begin, // position of the region with respect to the full molecule sequence
                    author_residue_number: m.str.begin, // position with respect to the beginning of the structure (in this case first three residues are not observed, i.e. residues 27, 28, 29 with respect to the full molecule)
                    author_insertion_code: undefined,
                },
                end: {
                    residue_number: m.seq.end,
                    author_residue_number: m.str.end,
                    author_insertion_code: undefined,
                }
            })
        }

        return mappingJson;
    });
}

function filterAnnotations(tabOrig) {
    const tab = _.cloneDeep(tabOrig);
    console.log('tab', tab);

    Object.values(annotationFields).filter(af => af.variant).map(af => af.name).forEach(annotationName => {
        tab[annotationName] = tab[annotationName].map(val => {
            if (val === "0" || val === "no") return undefined;
            if (val === "yes") return 1;
            return val;
        })
    })

    // [annotationFields.gscount.name, annotationFields.pscount.name].forEach(annotationName => {
    //     tab[annotationName] = tab[annotationName].map(val => val === "0" ? undefined : val);
    // })

    return tab;
}

function createDerivedAnnotations(tabOrig) {
    const tab = _.cloneDeep(tabOrig);
    let per_3d = [];
    for (let i = 0; i < tab['pvalue'].length; ++i ){
        per_3d.push(tab['OR'][i] > 1 && tab['pvalue'][i] < 0.05);
    }
    tab['PER_3D'] = per_3d;

    return tab;
}

function getFeatures(annotationsFileName, filters) {

    function createColorScale(af, vals){

        let colorScale;
        let colors;
        let distinctVals;
        const scale = chroma.scale([chroma(af.color).brighten(3).hex(), af.color]);
        if (af.categorical) {
            distinctVals = Array.from(new Set(vals))
            colors = scale.colors(distinctVals);                        
            return val => colors[distinctVals.indexOf(val)];;

        } else {
            const valsFiltered = vals.filter(val => !(isNaN(parseInt(val))))
            
            if (valsFiltered.length === 0) return () => chroma("white").hex();

            const min  = Math.min(...valsFiltered);
            const max  = Math.max(...valsFiltered);

            colorScale = scale.domain([min,max]);

            return val => isNaN(parseInt(val)) ? chroma("white").hex() : colorScale(val).hex();
        }
        

        

        // if (annotation.name === 'PVE3D') {
        //     //return c => c ? "#FF0000": "#FFFFFF";
        //     return c => annotation.color;
        // } else if (annotationName === annotationFields.pscount.name) {
        //     return c => "#FF0000";
        // } else if (annotationName === annotationFields.gscount.name) {
        //     return c => "#0000FF";
        // } else {
        //     let colorScale;
        //     let colors;
        //     let distinctVals;
        //     if (af.categorical) {
        //         distinctVals = Array.from(new Set(vals));
        //         if (annotationName === annotationFields.HotSpot3D.name){
        //             colors = chroma.scale(['gray', 'red']).colors(2);
        //         } else {
        //             colors = chroma.scale('RdYlBu').colors(distinctVals.length);
        //         }

        //     } else {
        //         const valsFiltered = vals.filter(val => !(isNaN(parseInt(val))))
        //         const min  = Math.min(...valsFiltered);
        //         const max  = Math.max(...valsFiltered);
        //         if (annotationName === annotationFields.paraz3dscore.name){
        //             colorScale = chroma.scale(['white', 'magenta']).domain([0,max]);
        //         } else if (annotationName === annotationFields.mtr3dscore.name){
        //             colorScale = chroma.scale(['purple', 'white']).domain([0,max]);
        //         } else {
        //             colorScale = chroma.scale(['#f00', '#0f0']).domain([min,max]);
        //         }
        //     }

        //     return function (val) {

        //         if (categorical) {
        //             return colors[distinctVals.indexOf(val)];

        //         } else {

        //             return isNaN(parseInt(val)) ? chroma("black").hex() : colorScale(val).hex() ;
        //         }
        //     }

        // }


    }
    return $.get(annotationsFileName).then(data => {
        let tab = parseTSV(data, filters);
        tab = filterAnnotations(tab);
        //tab = createDerivedAnnotations(tab); //PER_3D is in v2 explicitely stored in the tables as the PVE3D column

        const features = [];
        //Object.keys(annotationFields)
        Object.values(annotationFields)
            // .filter(afk => !annotationFields[afk].derived)
            .forEach(af => {
                const key = af.name;
                const vals = tab[key];
                const catName = af.category;

                let colorScale = createColorScale(af, vals);

                let lastIx = 0;
                let lastVal = tab[key][0];
                //merge same value neighboring annotations
                tab[key].forEach((val, ix) => {
                    if (ix !== 0 && lastVal !== val) {
                        // output only values which are defined
                        if (lastVal !== undefined) {
                            features.push({
                                type: key,
                                category: catName,
                                begin: String(tab['Uniprot_position'][lastIx]),
                                end: String(tab['Uniprot_position'][ix - 1]),
                                color: colorScale(lastVal),
                                // description: `${annotationFields[key].label}: ${lastVal}`
                                description: `${lastVal}`
                            });
                        }

                        lastIx = ix;
                        lastVal = val;
                    }
                });

                const val = tab[key][tab[key].length-1]
                if (val !== undefined) {
                    features.push({
                        type: key,
                        category: catName,
                        begin: String(tab['Uniprot_position'][lastIx]),
                        end: String(tab['Uniprot_position'][tab[key].length - 1]),
                        color: colorScale(val),
                        // description: `${annotationFields[key].label}: ${val}`
                        description: `${val}`
                    });
                }
        });

        return {
            sequence: tab['aminoacid'].join(''),
            features: features,
            data: tab
        };
    });
}

function getDataUris(params) {

    const uris = new DataUris();

    if (params.method === METHOD.PDB.id) {
        uris.annotations[params.structureId] = `${dataDir}pdb/Experimentally_solved_sinlge_${params.geneName}_${params.structureId}.txt`;
    } else if (params.method === METHOD.PDB_MULTI.id) {
        uris.annotations[params.structureId] = `${dataDir}pdb-multi/Experimentally_solved_complex_${params.structureId.split("_")[0]}.txt`;    
    } else if (params.method === METHOD.ALPHAFOLD.id) {
        uris.annotations[params.structureId] = `${dataDir}/alphafold/Alpha_fold_${params.geneName}.txt`;
        uris.structure = `${dataDir}alphafold/${params.geneName}.pdb`;
    }

    return uris;
}

function generateCustomConfig() {

    let customConfig = {
        "categories": [
            {
                name: "NDD",
                label: "NDD",
                visualizationType: "basic"
            },
            {
                "name": "ANTIGEN",
                "label": "Antigen",
                "visualizationType": "basic"
            },
            {
                "name": "DOMAINS_AND_SITES",
                "label": "Domains & sites",
                "visualizationType": "basic"
            },
            {
                "name": "MOLECULE_PROCESSING",
                "label": "Molecule processing",
                "visualizationType": "basic"
            },
            {
                "name": "PTM",
                "label": "PTM",
                "visualizationType": "basic"
            },
            {
                "name": "SEQUENCE_INFORMATION",
                "label": "Sequence information",
                "visualizationType": "basic"
            },
            {
                "name": "STRUCTURAL",
                "label": "Structural features",
                "visualizationType": "basic"
            },
            {
                "name": "TOPOLOGY",
                "label": "Topology",
                "visualizationType": "basic"
            },
            {
                "name": "MUTAGENESIS",
                "label": "Mutagenesis",
                "visualizationType": "basic"
            },
            {
                "name": "PROTEOMICS",
                "label": "Proteomics",
                "visualizationType": "basic"
            }
            ,{
                "name": "VARIATION",
                "label": "Variants",
                "visualizationType": "variant"
            }
        ],
        "trackNames": {
            "chain": {
                "label": "Chain",
                "tooltip": "(aka mature region). This describes the extent of a polypeptide chain in the mature protein following processing"
            },
            "transit": {
                "label": "Transit peptide",
                "tooltip": "This describes the extent of a transit peptide"
            },
            "init_met": {
                "label": "Initiator methionine",
                "tooltip": "This indicates that the initiator methionine is cleaved from the mature protein"
            },
            "propep": {
                "label": "Propeptide",
                "tooltip": "Part of a protein that is cleaved during maturation or activation"
            },
            "peptide": {
                "label": "Peptide",
                "tooltip": "The position and length of an active peptide in the mature protein"
            },
            "signal": {
                "label": "Signal",
                "tooltip": "N-terminal signal peptide"
            },
            "helix": {
                "label": "Helix",
                "tooltip": "The positions of experimentally determined helical regions"
            },
            "strand": {
                "label": "Beta strand",
                "tooltip": "The positions of experimentally determined beta strands"
            },
            "turn": {
                "label": "Turn",
                "tooltip": "The positions of experimentally determined hydrogen-bonded turns"
            },
            "disulfid": {
                "label": "Disulfide bond",
                "tooltip": "The positions of cysteine residues participating in disulphide bonds"
            },
            "crosslnk": {
                "label": "Cross-link",
                "tooltip": "Covalent linkages of various types formed between two proteins or between two parts of the same protein"
            },
            "region": {
                "label": "Region",
                "tooltip": "Regions in multifunctional enzymes or fusion proteins, or characteristics of a region, e.g., protein-protein interactions mediation"
            },
            "coiled": {
                "label": "Coiled coil",
                "tooltip": "Coiled coils are built by two or more alpha-helices that wind around each other to form a supercoil"
            },
            "motif": {
                "label": "Motif",
                "tooltip": "Short conserved sequence motif of biological significance"
            },
            "repeat": {
                "label": "Repeat",
                "tooltip": "Repeated sequence motifs or repeated domains within the protein"
            },
            "ca_bind": {
                "label": "Calcium binding",
                "tooltip": "Calcium-binding regions, such as the EF-hand motif"
            },
            "dna_bind": {
                "label": "DNA binding",
                "tooltip": "DNA-binding domains such as AP2/ERF domain, the ETS domain, the Fork-Head domain, the HMG box and the Myb domain"
            },
            "domain": {
                "label": "Domain",
                "tooltip": "Specific combination of secondary structures organized into a characteristic three-dimensional structure or fold"
            },
            "zn_fing": {
                "label": "Zinc finger",
                "tooltip": "Small, functional, independently folded domain that coordinates one or more zinc ions"
            },
            "np_bind": {
                "label": "Nucleotide binding",
                "tooltip": "(aka flavin-binding). Region in the protein which binds nucleotide phosphates"
            },
            "metal": {
                "label": "Metal binding",
                "tooltip": "Binding site for a metal ion"
            },
            "site": {
                "label": "Site",
                "tooltip": "Any interesting single amino acid site on the sequence"
            },
            "binding": {
                "label": "Binding site",
                "tooltip": "Binding site for any chemical group (co-enzyme, prosthetic group, etc.)"
            },
            "act_site": {
                "label": "Active site",
                "tooltip": "Amino acid(s) directly involved in the activity of an enzyme"
            },
            "mod_res": {
                "label": "Modified residue",
                "tooltip": "Modified residues such as phosphorylation, acetylation, acylation, methylation"
            },
            "lipid": {
                "label": "Lipidationasdf",
                "tooltip": "Covalently attached lipid group(s)"
            },
            "carbohyd": {
                "label": "Glycosylation",
                "tooltip": "Covalently attached glycan group(s)"
            },
            "compbias": {
                "label": "Compositional bias",
                "tooltip": "Position of regions of compositional bias within the protein and the particular amino acids that are over-represented within those regions"
            },
            "conflict": {
                "label": "Sequence conflict",
                "tooltip": "Sequence discrepancies of unknown origin"
            },
            "non_cons": {
                "label": "Non-adjacent residues",
                "tooltip": "Indicates that two residues in a sequence are not consecutive and that there is an undetermined number of unsequenced residues between them"
            },
            "non_ter": {
                "label": "Non-terminal residue",
                "tooltip": "The sequence is incomplete. The residue is not the terminal residue of the complete protein"
            },
            "unsure": {
                "label": "Sequence uncertainty",
                "tooltip": "Regions of a sequence for which the authors are unsure about the sequence assignment"
            },
            "non_std": {
                "label": "Non-standard residue",
                "tooltip": "Non-standard amino acids (selenocysteine and pyrrolysine)"
            },
            "mutagen": {
                "label": "Mutagenesis",
                "tooltip": "Site which has been experimentally altered by mutagenesis"
            },
            "topo_dom": {
                "label": "Topological domain",
                "tooltip": "Location of non-membrane regions of membrane-spanning proteins"
            },
            "transmem": {
                "label": "Transmembrane",
                "tooltip": "Extent of a membrane-spanning region"
            },
            "intramem": {
                "label": "Intramembrane",
                "tooltip": "Extent of a region located in a membrane without crossing it"
            },
            "variant": {
                "label": "Natural variant",
                "tooltip": "Natural variant of the protein, including polymorphisms, variations between strains, isolates or cultivars, disease-associated mutations and RNA editing events"
            },
            "unique": {
                "label": "Unique peptide",
                "tooltip": ""
            },
            "non_unique": {
                "label": "Non-unique peptide",
                "tooltip": ""
            }
        }
    };

    Object.values(annotationFields).forEach(af => {
        customConfig.trackNames[af.name.toLocaleLowerCase()] = {
            label: af.label,
            tooltip:af.tooltip
        }
    });

    return 'data:text/plain;charset=utf-8,' + encodeURI(JSON.stringify(customConfig));

}

function convertToLiteMolRgb(chromaRgb){
    return {
        r: chromaRgb[0]/255.0,
        g: chromaRgb[1]/255.0,
        b: chromaRgb[2]/255.0,
    }
}

function getExtraHighlights(tab) {

    const variantsContent = [];
    Object.values(annotationFields).filter(af => af.variant).forEach(af => {
        const variants = tab[af.name].map((v, i) => v === undefined ? undefined : parseInt(tab.Uniprot_position[i])).filter(v => v !== undefined);
        variantsContent.push(
            {
                sequenceNumbers: variants,
                atomNames: ['CA'],
                label: af.label,
                visual: {
                    type: 'BallsAndSticks',
                    params: {useVDW: true, vdwScaling: 1.2, bondRadius: 0.13, detail: 'Automatic'},
                    color: convertToLiteMolRgb(chroma(af.color).rgb()),
                    alpha: 1
                }
    
            }
        )

    })

    var extraHighlights = {
        controlVisibility: true, //whether the list of custom highlights will be shown as a dropdown through which the can control visibility of the individual highlights
        label: 'Variation highlighting',
        content: variantsContent
    }

    return extraHighlights;
}

function startMolArtWithFeatures(params) {

    // const uniprotId = sequenceSelection.select2('data')[0].id;

    return getFeatures(params.dataUris.annotations[params.structureId], params.filters).then(features => {
        const molArtParams = {
            uniprotId: params.uniprotId,
            containerId: containerId,
            alwaysLoadPredicted: false,
            customDataSources : [{
                source: 'ES-NDD',
                useExtension: false,
                data: features
            }],
            exclusions: ['PREDICT_PROTEIN', 'VARIATION', 'ANTIGEN', 'MUTAGENESIS', 'PROTEOMICS'],
            extraHighlights: getExtraHighlights(features.data),
            lmInitSurfaceTransparency: 25,
            customConfig: customConfig
        };

        if (params.ssMapping) {
            molArtParams.sequenceStructureMapping = [params.ssMapping];
        }

        if (params.pdbIds) {
            molArtParams.pdbIds = params.pdbIds;
        }
        window.molartPlugin = new MolArt(molArtParams);
        window.molartPlugin.on("pvReady", () => {
            $('.up_pftv_category-name[title="USER_PROVIDED_STRUCTURES"]').text("Selected structure");
            $('.up_pftv_category-name[title="EXPERIMENTAL_STRUCTURES"]').text("Selected structure");
        })
    })

}

const btnShowOnClick = function(data) {
    if ($('#btnShow').hasClass("disabled")) {
        // $('#warning').css("display", "block");
        return;
    }

    const gene = geneSelection.select2('data')[0].id;
    const method = methodSelection.select2('data')[0].id;
    const structureId = structureSelection.select2('data')[0].id;

    const uniprotId = Object.keys(esNddData[gene][method])[0];

    const $container = $(`#${containerId}`);
    $container.empty();
    if (window.molartPlugin) {
        try {
            window.molartPlugin.destroy();
        } catch (e) {
        }
    }

    let molArtPromise;
    if (method === METHOD.PDB.id){
        // const pdbIds = data[gene][method][uniprotId].map(pdbId => pdbId.replace('_', ':'));
        const pdbIds = [structureId.replace('_', ':')];
        const dataUris = getDataUris({
            geneName: gene,
            method: method,
            structureId: structureId
        });
        molArtPromise = startMolArtWithFeatures({uniprotId: uniprotId, dataUris: dataUris, structureId: structureId, pdbIds:pdbIds});
    }else if (method === METHOD.PDB_MULTI.id){
        // const pdbIds = data[gene][method][uniprotId].map(pdbId => pdbId.replace('_', ':'));
        const sStructureId = structureId.split("_");
        const pdbIds = [structureId.replace('_', ':')];
        const dataUris = getDataUris({
            geneName: gene,
            method: method,
            structureId: structureId
        });
        molArtPromise = startMolArtWithFeatures({
            uniprotId: uniprotId,
            dataUris: dataUris,
            structureId: structureId,
            pdbIds:pdbIds,
            filters:{
                gene: gene,
                chain: sStructureId[1]
            }});
    } else if (method === METHOD.ALPHAFOLD.id){

        const dataUris = getDataUris({
            geneName: gene,
            method: method,
            structureId: structureId
        });

        let chainId = 'A';
        molArtPromise = getSSMapping({
            mappingUri: dataUris.annotations[structureId],
            structureUri: dataUris.structure,
            pdbId: structureId,
            chainId: chainId
        }).then(mappingsJson => {
            startMolArtWithFeatures({
                uniprotId: uniprotId,
                dataUris: dataUris,
                structureId: structureId,
                ssMapping: mappingsJson
            });
        })

    } 

    return molArtPromise.then(() => {
        $("#molartWrapperContainer").css("display", "block" );
        $("#mainInfoContainer").css("display", "none");
        $([document.documentElement, document.body]).animate({
            scrollTop: $container.offset().top
        }, 1000);
    });

    // const pdbUri = getPdbUri(gene, method, sequence);

};

// function loadData() {
//     return $.getJSON('data/list.json')
// }

function populateGeneList(data) {

    function sequenceChange(sequence) {
        structureSelection.html("");
        const gene = geneSelection.select2('data')[0].id;
        const method = methodSelection.select2('data')[0].id;
        structureSelection.val(null).trigger('change');
        data[gene][method][sequence].sort().forEach((structure, ix) => {
            structureSelection.append(new Option(structure, structure, ix === 0, ix === 0));
        });
        $('#btnShow').removeClass("disabled");
    }

    sequenceSelection.on('select2:select', function (e) {
        sequenceChange(e.params.data.id);
    });

    function methodChange(method) {
        sequenceSelection.html("");
        const gene = geneSelection.select2('data')[0].id;
        sequenceSelection.val(null).trigger('change');
        let firstSequenceId;
        Object.keys(data[gene][method]).forEach((sequence, ix) => {
            sequenceSelection.append(new Option(sequence, sequence, ix === 0, ix === 0));
            if (ix === 0) firstSequenceId = sequence;
        });
        sequenceChange(firstSequenceId);
    }

    methodSelection.on('select2:select', function (e) {
        methodChange(e.params.data.id);
    });

    geneSelection.on('select2:select', function (e) {
        const gene = e.params.data.id;
        if (Object.keys(data).indexOf(gene) === -1) return;

        methodSelection.html("");
        let firstMethodId;
        Object.keys(data[gene]).sort().forEach((method, ix)   => {
            const opt = methodSelection.append(new Option(METHOD[method].label, METHOD[method].id, ix === 0, ix === 0));
            if (ix === 0) firstMethodId = method;
        });
        methodChange(firstMethodId);

        // let firstSequenceId;
        // Object.keys(data[gene][firstMethodId]).forEach((sequence, ix)   => {
        //     const opt = sequenceSelection.append(new Option(sequence, sequence, ix === 0, ix === 0));
        //     if (ix === 0) firstSequenceId = sequence;
        // });
        // sequenceChange(firstSequenceId);
    });

    Object.keys(data).sort().forEach((gene) => {
        var newOption = new Option(gene, gene, false, false);
        geneSelection.append(newOption);
    });
}