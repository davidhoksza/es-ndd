
import 'bootstrap/dist/js/bootstrap.bundle.min'
import _ from 'lodash'
import chroma from 'chroma-js'
import $ from 'jquery'
import 'select2'
import MolArt from 'molart'

window.$ = $

import 'bootstrap/dist/css/bootstrap.min.css'
import '@fortawesome/fontawesome-free/css/all.css'
import 'simple-line-icons/css/simple-line-icons.css'
import 'select2/dist/css/select2.min.css'
import '../css/ndd.css'

import esNddData from '../../data/list.json'

const METHOD = {
    PDB: {
        id: 'PDB',
        label: 'Monomeric structure from Protein Data Bank'
    },
    PDB_MULTI: {
        id: 'PDB_MULTI',
        label: 'Complex structure from Protein Data Bank'
    },
    SWISS:
        {
            id: 'SWISS',
            label: 'Homology models from SWISS-PROT'
        },
    RAPTOR: {
        id: "RAPTOR",
        label: 'Predicted structures from RaptorX'
    }
};

const containerId = 'molartContainer';

const annotationFields = {
    HotSpot3D: {
        name: "HotSpot3D",
        label: 'Essential sites',
        tooltip: 'Essential sites',
        category: "Structural based annotations"
        , categorical: true
        , derived: false
    },
    paraz3dscore: {
        name: "paraz3dscore",
        label: 'Paralog conserved sites',
        tooltip: 'Paralog conserved sites',
        category: "Structural based annotations"
        , categorical: false
        , derived: false
    },
    mtr3dscore: {
        name: "mtr3dscore",
        label: 'Missense constraint sites',
        tooltip: 'Missense constraint sites',
        category: "Structural based annotations"
        , categorical: false
        , derived: false
    },
    PER_3D: {
        name: "PER_3D",
        label: 'Variant enriched sites',
        tooltip: 'Variant enriched sites',
        category: "Structural based annotations"
        , categorical: true
        , derived: true
    },
    gscount: {
        name: "gscount",
        label: 'Neutral variants',
        tooltip: 'Neutral variants',
        category: "Variants"
        , categorical: false
        , derived: false
    },
    pscount: {
        name: "pscount",
        label: 'Pathogenic variants',
        tooltip: 'Pathogenic variants',
        category: "Variants"
        , categorical: false
        , derived: false
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
    sequenceSelection.select2({
        placeholder: "Select a sequence",
        allowClear: false
    });
    structureSelection.select2({
        placeholder: "Select a structure",
        allowClear: false
    });
    // return loadData().then(data => {
    //     console.log('data', data);
    //     populateGeneList(data);
    //     $('#btnShow').on("click", () => btnShowOnClick(data))
    // });
    console.log('data', esNddData);
    populateGeneList(esNddData);
    $('#btnShow').on("click", () => btnShowOnClick(esNddData))
});

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

    return tab;
}

// Get sequence-structure mapping
function getSSMapping(params) {
    return $.get(params.mappingUri).then(data => {
        const tab = parseTSV(data, params.filters);
        const seqIxs = tab['Uniprot_position'].map(d => parseInt(d));
        const strIxs = tab['Position_in_structure'].map(d => parseInt(d));

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

    [annotationFields.gscount.name, annotationFields.pscount.name].forEach(annotationName => {
        tab[annotationName] = tab[annotationName].map(val => val === "0" ? undefined : val);
    })

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

    function createColorScale(annotationName, vals, categorical){

        if (annotationName === 'PER_3D') {
            return c => c ? "#FF0000": "#FFFFFF";
        } else if (annotationName === annotationFields.pscount.name) {
            return c => "#FF0000";
        } else if (annotationName === annotationFields.gscount.name) {
            return c => "#0000FF";
        } else {
            let colorScale;
            let colors;
            let distinctVals;
            if (categorical) {
                distinctVals = Array.from(new Set(vals));
                if (annotationName === annotationFields.HotSpot3D.name){
                    colors = chroma.scale(['gray', 'red']).colors(2);
                } else {
                    colors = chroma.scale('RdYlBu').colors(distinctVals.length);
                }

            } else {
                const valsFiltered = vals.filter(val => !(isNaN(parseInt(val))))
                const min  = Math.min(...valsFiltered);
                const max  = Math.max(...valsFiltered);
                if (annotationName === annotationFields.paraz3dscore.name){
                    colorScale = chroma.scale(['white', 'magenta']).domain([0,max]);
                } else if (annotationName === annotationFields.mtr3dscore.name){
                    colorScale = chroma.scale(['purple', 'white']).domain([0,max]);
                } else {
                    colorScale = chroma.scale(['#f00', '#0f0']).domain([min,max]);
                }
            }

            return function (val) {

                if (categorical) {
                    return colors[distinctVals.indexOf(val)];

                } else {

                    return isNaN(parseInt(val)) ? chroma("black").hex() : colorScale(val).hex() ;
                }
            }

        }


    }
    return $.get(annotationsFileName).then(data => {
        let tab = parseTSV(data, filters);
        tab = filterAnnotations(tab);
        tab = createDerivedAnnotations(tab);

        const features = [];
        Object.keys(annotationFields)
            // .filter(afk => !annotationFields[afk].derived)
            .forEach(key => {
            const vals = tab[key];
            const catName = annotationFields[key].category;

            let colorScale = createColorScale(key, vals, annotationFields[key].categorical);

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
                            description: `${annotationFields[key].label}: ${lastVal}`
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
                    description: `${annotationFields[key].label}: ${val}`
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
        uris.annotations[params.structureId] = `data/pdb/Experimentally_solved_sinlge_${params.geneName}_${params.structureId}.txt`;
    } else if (params.method === METHOD.PDB_MULTI.id) {
        uris.annotations[params.structureId] = `data/pdb-multi/Experimentally_solved_complex_${params.structureId.split("_")[0]}.txt`;
    } else if (params.method === METHOD.SWISS.id) {
        uris.annotations[params.structureId] = `data/swissmodel/Swiss_${params.structureId}.txt`;
        uris.structure = `data/swissmodel/${params.geneName}.pdb`;
    } else if (params.method === METHOD.RAPTOR.id) {
        uris.annotations[params.structureId] = `data/raptor/RaptorX_${params.geneName}.txt`;
        uris.structure = `data/raptor/${params.geneName}.pdb`;
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

    Object.keys(annotationFields).forEach(k => {
        customConfig.trackNames[k.toLocaleLowerCase()] = {
            label: annotationFields[k].label,
            tooltip:  annotationFields[k].tooltip
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
    const patientVariants = tab[annotationFields.pscount.name].map((v, i) => v === undefined ? undefined : tab.Position_in_structure[i]).filter(v => v !== undefined);
    const populationVariants = tab[annotationFields.gscount.name].map((v, i) => v === undefined ? undefined : tab.Position_in_structure[i]).filter(v => v !== undefined);

    var extraHighlights = {
        controlVisibility: true, //whether the list of custom highlights will be shown as a dropdown through which the can control visibility of the individual highlights
        label: 'Variation highlighting',
        content: [{
            sequenceNumbers: patientVariants,
            atomNames: ['CA'],
            label: 'Pathogenic',
            visual: {
                type: 'BallsAndSticks',
                params: {useVDW: true, vdwScaling: 1.2, bondRadius: 0.13, detail: 'Automatic'},
                color: convertToLiteMolRgb(chroma("red").rgb()),
                alpha: 1
            }

        }, {
            sequenceNumbers: populationVariants,
            atomNames: ['CA'],
            label: 'Neutral',
            visual: {
                type: 'BallsAndSticks',
                params: {useVDW: true, vdwScaling: 1.2, bondRadius: 0.13, detail: 'Automatic'},
                color: convertToLiteMolRgb(chroma("blue").rgb()),
                alpha: 1
            }
        }, {
            sequenceNumbers: _.intersection(populationVariants, patientVariants),
            atomNames: ['CA'],
            label: 'Both',
            visual: {
                type: 'BallsAndSticks',
                params: {useVDW: true, vdwScaling: 1.2, bondRadius: 0.13, detail: 'Automatic'},
                color: convertToLiteMolRgb(chroma("orange").rgb()),
                alpha: 1
            }
        }]
    }

    return extraHighlights;
}

function startMolArtWithFeatures(params) {

    const uniprotId = sequenceSelection.select2('data')[0].id;

    return getFeatures(params.dataUris.annotations[params.structureId], params.filters).then(features => {
        const molArtParams = {
            uniprotId: uniprotId,
            containerId: containerId,
            alwaysLoadPredicted: false,
            customDataSources : [{
                source: 'ES-NDD',
                useExtension: false,
                data: features
            }],
            exclusions: ['PREDICT_PROTEIN', 'VARIATION', 'ANTIGEN', 'MUTAGENESIS', 'PROTEOMICS'],
            extraHighlights: getExtraHighlights(features.data),
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
        molArtPromise = startMolArtWithFeatures({dataUris: dataUris, structureId: structureId, pdbIds:pdbIds});
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
            dataUris: dataUris,
            structureId: structureId,
            pdbIds:pdbIds,
            filters:{
                gene: gene,
                chain: sStructureId[1]
            }});
    } else if (method === METHOD.RAPTOR.id){

        const dataUris = getDataUris({
            geneName: gene,
            method: method,
            structureId: structureId
        });

        let chainId = ' ';
        molArtPromise = getSSMapping({
            mappingUri: dataUris.annotations[structureId],
            structureUri: dataUris.structure,
            pdbId: structureId,
            chainId: chainId
        }).then(mappingsJson => {
            startMolArtWithFeatures({dataUris: dataUris, structureId: structureId, ssMapping:mappingsJson});
        })

    } else if (method === METHOD.SWISS.id){

        const sStructureId = structureId.split("_");
        const pdbId = sStructureId[0];
        const chainId = sStructureId[1];

        const dataUris = getDataUris({
            geneName: gene,
            structureId: structureId,
            chainId: chainId,
            method: method,
        });

        molArtPromise = getSSMapping({
            mappingUri: dataUris.annotations[structureId],
            structureUri: dataUris.structure,
            pdbId:pdbId,
            chainId: chainId
        }).then(mappingsJson => {
            startMolArtWithFeatures({dataUris: dataUris, structureId: structureId, ssMapping:mappingsJson});
        });
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
        data[gene][method][sequence].forEach((structure, ix) => {
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

        let firstSequenceId;
        Object.keys(data[gene][firstMethodId]).forEach((sequence, ix)   => {
            const opt = sequenceSelection.append(new Option(sequence, sequence, ix === 0, ix === 0));
            if (ix === 0) firstSequenceId = sequence;
        });
        sequenceChange(firstSequenceId);
    });

    Object.keys(data).sort().forEach((gene) => {
        var newOption = new Option(gene, gene, false, false);
        geneSelection.append(newOption);
    });
}