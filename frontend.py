import json
import os
from io import StringIO

import pandas as pd
import plotly.express as px
import requests
import streamlit as st

API_URL = os.getenv("API_URL", "http://localhost:8000")

st.set_page_config(
    page_title="Microbiome Analysis Dashboard",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="expanded",
)


st.markdown(
    """
<style>
    .main-header {
        font-size: 2.5rem;
        color: #4B89DC;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #656565;
        margin-bottom: 0.5rem;
    }
    .success-message {
        color: #28A745;
        font-weight: bold;
        padding: 0.5rem;
        border-radius: 0.3rem;
        background-color: #E6F4EA;
    }
    .info-text {
        color: #656565;
        font-size: 0.9rem;
    }
</style>
""",
    unsafe_allow_html=True,
)


st.markdown(
    "<h1 class='main-header'>Microbiome Analysis Dashboard</h1>", unsafe_allow_html=True
)


def upload_file(file, file_type):
    """Upload file to API and return the parsed data"""
    if file_type == "csv":
        files = {"file": (file.name, file.getvalue(), f"application/{file_type}")}
        response = requests.post(f"{API_URL}/upload/", files=files)
    elif file_type == "json":
        try:

            json_data = json.loads(file.getvalue().decode("utf-8"))

            if (
                "feature_ids" not in json_data or "counts_matrix" not in json_data
            ) and ("taxonomy_ids" in json_data or "abundance_matrix" in json_data):

                transformed_data = {
                    "sample_ids": json_data.get("sample_ids", []),
                    "feature_ids": json_data.get(
                        "taxonomy_ids", json_data.get("feature_ids", [])
                    ),
                    "counts_matrix": json_data.get(
                        "abundance_matrix", json_data.get("counts_matrix", [])
                    ),
                }

                if "metadata" in json_data:
                    transformed_data["metadata"] = json_data["metadata"]

                response = requests.post(
                    f"{API_URL}/upload/",
                    files={
                        "json_data": (
                            None,
                            json.dumps(transformed_data),
                            "application/json",
                        )
                    },
                )
            else:
                response = requests.post(
                    f"{API_URL}/upload/",
                    files={
                        "json_data": (None, json.dumps(json_data), "application/json")
                    },
                )
        except json.JSONDecodeError as e:
            st.error(f"Error parsing JSON: {str(e)}")
            return None
        except Exception as e:
            st.error(f"Error processing JSON file: {str(e)}")
            return None
    else:
        st.error("Unsupported file type")
        return None

    if response.status_code == 200:
        return response.json()
    else:
        st.error(f"Error uploading file: {response.text}")
        return None


def perform_diversity_analysis(data, metrics):
    """Send data to diversity analysis endpoint"""

    body = {"request": {"metrics": metrics}, "data": data}

    response = requests.post(
        f"{API_URL}/diversity/analyze",
        json=body,
        headers={"Content-Type": "application/json"},
    )

    if response.status_code == 200:
        return response.json()
    else:
        st.error(f"Error in diversity analysis: {response.text}")
        return None


def perform_correlation_analysis(
    data, correlation_method, metadata_columns, min_abundance, min_prevalence
):
    """Send data to correlation analysis endpoint"""

    upload_response = requests.post(
        f"{API_URL}/upload/",
        files={"json_data": (None, json.dumps(data), "application/json")},
    )

    if upload_response.status_code != 200:
        st.error(f"Error setting current data: {upload_response.text}")
        return None

    cleaned_metadata_columns = None
    if metadata_columns:
        cleaned_metadata_columns = [
            col.replace("metadata_", "") for col in metadata_columns
        ]

    request_data = {
        "correlation_method": correlation_method,
        "metadata_columns": cleaned_metadata_columns,
        "min_abundance": min_abundance,
        "min_prevalence": min_prevalence,
    }

    response = requests.post(
        f"{API_URL}/correlation/analyze",
        json=request_data,
        headers={"Content-Type": "application/json"},
    )

    if response.status_code == 200:
        return response.json()
    else:
        st.error(f"Error in correlation analysis: {response.text}")
        return None


def plot_alpha_diversity(alpha_diversity, sample_ids):
    """Create alpha diversity visualizations"""
    st.markdown(
        "<h2 class='sub-header'>Alpha Diversity Metrics</h2>", unsafe_allow_html=True
    )

    tab1, tab2 = st.tabs(["Bar Charts", "Heatmap"])

    with tab1:
        for metric, values in alpha_diversity.items():
            metric_values = [values[sample_id] for sample_id in sample_ids]
            df = pd.DataFrame(
                {"Sample": sample_ids, f"{metric.capitalize()} Index": metric_values}
            )

            fig = px.bar(
                df,
                x="Sample",
                y=f"{metric.capitalize()} Index",
                title=f"{metric.capitalize()} Diversity Index by Sample",
                color=f"{metric.capitalize()} Index",
                color_continuous_scale="Viridis",
            )
            fig.update_layout(height=500)
            st.plotly_chart(fig, use_container_width=True)

    with tab2:

        heatmap_data = pd.DataFrame(index=sample_ids)

        for metric, values in alpha_diversity.items():
            heatmap_data[metric.capitalize()] = [
                values[sample_id] for sample_id in sample_ids
            ]

        normalized_data = (heatmap_data - heatmap_data.min()) / (
            heatmap_data.max() - heatmap_data.min()
        )

        fig = px.imshow(
            normalized_data.T,
            title="Alpha Diversity Metrics Heatmap (Normalized)",
            color_continuous_scale="Viridis",
            labels=dict(x="Sample", y="Metric", color="Normalized Value"),
        )
        fig.update_layout(height=600)
        st.plotly_chart(fig, use_container_width=True)


def plot_beta_diversity(beta_diversity, sample_ids):
    """Create beta diversity visualizations"""
    st.markdown(
        "<h2 class='sub-header'>Beta Diversity Metrics</h2>", unsafe_allow_html=True
    )

    for metric, distance_matrix in beta_diversity.items():
        st.markdown(f"#### {metric.capitalize()} Distance Matrix")

        df = pd.DataFrame(distance_matrix, index=sample_ids, columns=sample_ids)

        fig = px.imshow(
            df,
            title=f"{metric.capitalize()} Distance Matrix",
            labels=dict(x="Sample", y="Sample", color="Distance"),
            color_continuous_scale="Inferno",
        )
        fig.update_layout(height=600)
        st.plotly_chart(fig, use_container_width=True)

        if len(sample_ids) >= 3:
            from sklearn.manifold import MDS

            try:

                mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
                positions = mds.fit_transform(distance_matrix)

                mds_df = pd.DataFrame(
                    {
                        "Sample": sample_ids,
                        "Dimension 1": positions[:, 0],
                        "Dimension 2": positions[:, 1],
                    }
                )

                fig = px.scatter(
                    mds_df,
                    x="Dimension 1",
                    y="Dimension 2",
                    text="Sample",
                    title=f"{metric.capitalize()} MDS Plot",
                    labels={"Dimension 1": "Dimension 1", "Dimension 2": "Dimension 2"},
                )
                fig.update_traces(textposition="top center")
                fig.update_layout(height=600)
                st.plotly_chart(fig, use_container_width=True)
            except Exception as e:
                st.warning(f"Could not create MDS plot: {str(e)}")


def plot_correlation_matrix(correlation_data, title, threshold=0.5):
    """Create a correlation matrix heatmap"""

    taxa = list(correlation_data.keys())
    corr_matrix = pd.DataFrame(index=taxa, columns=taxa)

    for taxon1 in taxa:
        for taxon2 in taxa:
            if taxon2 in correlation_data[taxon1]:
                corr_matrix.loc[taxon1, taxon2] = correlation_data[taxon1][taxon2]
            else:
                corr_matrix.loc[taxon1, taxon2] = 0

    corr_matrix = corr_matrix.fillna(0)

    fig = px.imshow(
        corr_matrix,
        title=title,
        color_continuous_scale="RdBu_r",
        range_color=[-1, 1],
        labels=dict(x="Taxon", y="Taxon", color="Correlation"),
    )
    fig.update_layout(height=700, width=700)
    return fig


def plot_metadata_correlation(metadata_correlations, p_values, title):
    """Create a heatmap for metadata correlations"""
    taxa = list(metadata_correlations.keys())
    if not taxa:
        st.warning("No metadata correlations available")
        return None

    metadata_columns = set()
    for taxon_data in metadata_correlations.values():
        metadata_columns.update(taxon_data.keys())
    metadata_columns = list(metadata_columns)

    data = []
    for taxon in taxa:
        for meta_col in metadata_columns:
            if meta_col in metadata_correlations[taxon]:
                corr_value = metadata_correlations[taxon][meta_col]
                p_value = p_values[taxon].get(meta_col, 1.0) if p_values else 1.0
                data.append(
                    {
                        "Taxon": taxon,
                        "Metadata": meta_col,
                        "Correlation": corr_value,
                        "P-value": p_value,
                        "Significant": p_value < 0.05,
                    }
                )

    if not data:
        st.warning("No metadata correlations data available")
        return None

    df = pd.DataFrame(data)
    pivot_df = df.pivot(index="Taxon", columns="Metadata", values="Correlation")

    fig = px.imshow(
        pivot_df,
        title=title,
        color_continuous_scale="RdBu_r",
        color_continuous_midpoint=0,
        labels=dict(x="Metadata", y="Taxon", color="Correlation"),
    )
    fig.update_layout(height=600)

    return fig


with st.sidebar:
    st.title("Analysis Controls")

    st.header("Data Upload")
    uploaded_file = st.file_uploader(
        "Upload your microbiome data", type=["csv", "json"]
    )

    if uploaded_file is not None:
        file_type = uploaded_file.name.split(".")[-1].lower()

        st.header("Analysis Type")
        analysis_type = st.radio(
            "Select Analysis Type", ["Diversity Analysis", "Correlation Analysis"]
        )

        if analysis_type == "Diversity Analysis":
            st.header("Diversity Parameters")
            diversity_metrics = st.multiselect(
                "Select Diversity Metrics",
                ["shannon", "simpson", "pielou", "chao1"],
                default=["shannon", "simpson"],
            )

        elif analysis_type == "Correlation Analysis":
            st.header("Correlation Parameters")
            correlation_method = st.selectbox(
                "Correlation Method", ["spearman", "pearson", "kendall"], index=0
            )

            min_abundance = st.slider(
                "Minimum Abundance",
                min_value=0.0,
                max_value=0.2,
                value=0.01,
                step=0.01,
                help="Minimum relative abundance for a taxon to be included",
            )

            min_prevalence = st.slider(
                "Minimum Prevalence",
                min_value=0.0,
                max_value=1.0,
                value=0.1,
                step=0.05,
                help="Minimum fraction of samples where a taxon must be present",
            )

            if file_type == "csv":

                df = pd.read_csv(StringIO(uploaded_file.getvalue().decode("utf-8")))
                metadata_options = [
                    col for col in df.columns if col.startswith("metadata_")
                ]
                if metadata_options:
                    metadata_columns = st.multiselect(
                        "Metadata Columns for Correlation", metadata_options
                    )
                else:
                    metadata_columns = []
                    st.info("No metadata columns found in the CSV file.")
            else:
                metadata_columns = st.text_input(
                    "Metadata Columns (comma-separated)", ""
                )
                if metadata_columns:
                    metadata_columns = [
                        col.strip() for col in metadata_columns.split(",")
                    ]
                else:
                    metadata_columns = []

        run_analysis = st.button("Run Analysis")


if uploaded_file is not None:
    st.markdown(
        f"<p class='info-text'>Selected file: {uploaded_file.name}</p>",
        unsafe_allow_html=True,
    )

    file_type = uploaded_file.name.split(".")[-1].lower()

    with st.expander("Data Preview"):
        if file_type == "csv":
            df = pd.read_csv(StringIO(uploaded_file.getvalue().decode("utf-8")))
            st.dataframe(df.head(10))
        elif file_type == "json":
            try:
                json_data = json.loads(uploaded_file.getvalue().decode("utf-8"))
                st.json(json_data)
            except Exception:
                st.error("Could not parse JSON data for preview.")

    if "run_analysis" in locals() and run_analysis:
        st.markdown(
            "<div class='success-message'>Analysis started...</div>",
            unsafe_allow_html=True,
        )

        data = upload_file(uploaded_file, file_type)

        if data:

            if analysis_type == "Diversity Analysis":
                st.markdown(
                    "<h2 class='sub-header'>Diversity Analysis Results</h2>",
                    unsafe_allow_html=True,
                )

                with st.spinner("Running diversity analysis..."):
                    results = perform_diversity_analysis(data, diversity_metrics)

                    if results:
                        st.success("Diversity analysis completed successfully!")

                        plot_alpha_diversity(
                            results["alpha_diversity_metrics"], results["sample_ids"]
                        )

                        if results["beta_diversity_metrics"]:
                            plot_beta_diversity(
                                results["beta_diversity_metrics"], results["sample_ids"]
                            )

                        if results["group_comparison_metrics"] and any(
                            results["group_comparison_metrics"].values()
                        ):
                            st.markdown(
                                "<h2 class='sub-header'>Group Comparison Results</h2>",
                                unsafe_allow_html=True,
                            )

                            for metric, comparisons in results[
                                "group_comparison_metrics"
                            ].items():
                                if comparisons:
                                    st.markdown(
                                        f"#### {metric.capitalize()} Group Comparisons"
                                    )

                                    comparison_data = []
                                    for comp_name, comp_values in comparisons.items():
                                        row = {"Comparison": comp_name}
                                        row.update(comp_values)
                                        comparison_data.append(row)

                                    comp_df = pd.DataFrame(comparison_data)
                                    st.dataframe(comp_df)

                                    sig_comps = comp_df[comp_df["significant"]]
                                    if not sig_comps.empty:
                                        st.markdown("**Significant Comparisons:**")
                                        st.dataframe(sig_comps)

            elif analysis_type == "Correlation Analysis":
                st.markdown(
                    "<h2 class='sub-header'>Correlation Analysis Results</h2>",
                    unsafe_allow_html=True,
                )

                with st.spinner("Running correlation analysis..."):

                    meta_cols = metadata_columns if metadata_columns else None

                    results = perform_correlation_analysis(
                        data,
                        correlation_method,
                        meta_cols,
                        min_abundance,
                        min_prevalence,
                    )

                    if results:
                        st.success("Correlation analysis completed successfully!")

                        st.markdown("#### Filtered Taxa")
                        st.write(
                            f"Number of taxa after filtering: "
                            f"{len(results['filtered_taxa'])}"
                        )

                        st.markdown("#### Taxon Correlation Matrix")
                        corr_fig = plot_correlation_matrix(
                            results["taxon_correlations"],
                            f"Taxon Correlation Matrix ({correlation_method})",
                        )
                        st.plotly_chart(corr_fig, use_container_width=True)

                        if results.get("metadata_correlations"):
                            st.markdown("#### Metadata Correlations")
                            meta_fig = plot_metadata_correlation(
                                results["metadata_correlations"],
                                results["p_values"],
                                f"Metadata Correlation Matrix ({correlation_method})",
                            )
                            if meta_fig is not None:
                                st.plotly_chart(meta_fig, use_container_width=True)

                            st.markdown("#### Detailed Metadata Correlations")

                            meta_corr_list = []
                            for taxon in results["filtered_taxa"]:
                                if taxon in results["metadata_correlations"]:
                                    for meta_col, corr_val in results[
                                        "metadata_correlations"
                                    ][taxon].items():
                                        p_val = results["p_values"][taxon].get(
                                            meta_col, 1.0
                                        )
                                        meta_corr_list.append(
                                            {
                                                "Taxon": taxon,
                                                "Metadata": meta_col,
                                                "Correlation": corr_val,
                                                "P-value": p_val,
                                                "Significant": p_val < 0.05,
                                            }
                                        )

                            if meta_corr_list:
                                meta_corr_df = pd.DataFrame(meta_corr_list)
                                meta_corr_df = meta_corr_df.sort_values(
                                    by=["Significant", "P-value"],
                                    ascending=[False, True],
                                )
                                st.dataframe(meta_corr_df)

                                sig_corrs = meta_corr_df[meta_corr_df["Significant"]]
                                if not sig_corrs.empty:
                                    st.markdown(
                                        "**Significant Metadata Correlations:**"
                                    )
                                    st.dataframe(sig_corrs)

else:
    st.info("Please upload a CSV or JSON file containing microbiome data to begin.")

    with st.expander("File Format Instructions"):
        st.markdown(
            """
        ### CSV Format
        - First column: sample_id
        - Remaining columns: feature_id (e.g. taxa names)
        - Remaining cells: abundance counts
        - Optional metadata columns should have the prefix 'metadata_'
        ### JSON Format
        ```json
        {
            "sample_ids": ["sample1", "sample2", ...],
            "feature_ids": ["taxon1", "taxon2", ...],
            "counts_matrix": [
                [count1_1, count1_2, ...],
                [count2_1, count2_2, ...],
                ...
            ],
            "metadata": {
                "sample1": {"metadata_age": 30, "metadata_gender": "female", ...},
                "sample2": {"metadata_age": 25, "metadata_gender": "male", ...},
                ...
            }
        }
        ```
        """
        )


st.markdown("---")
st.markdown(
    "<p class='info-text'>Microbiome Data Analysis Dashboard</p>",
    unsafe_allow_html=True,
)
