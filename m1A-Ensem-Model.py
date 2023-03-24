import streamlit as st
from streamlit_option_menu import option_menu
from PIL import Image

selected2 = option_menu(None, ["Home", "Predictor", "Dataset", "Citations"],
                        icons=['house', 'search', 'list-task', 'book'],
                        menu_icon="cast", default_index=0, orientation="horizontal")
# selected2

if selected2 == "Home":
    image = Image.open('Methodology.png')

    st.header("m1A-Ensem: Accurate Identification of 1-Methyladnosine Sites through Ensemble Models")
    st.image(image)

    #st.subheader("A web-server for the prediction of 1-methyladenosine in RNA"
     #            "modifications.")
    st.write("1-methyladenosine (m1A) is a variant of methyladenosine that holds a methyl substituent in the 1st position. This modification has a prominent role in RNA stability and human metabolites. It is functionally related to adenosine. Traditional approaches, such as mass spectrometry and site-directed mutagenesis, proved to be time-consuming and complicated. The systematic screening of transformed sites using RNA sequences is gaining popularity at the moment. The present research focused on the identification of m1A sites within RNA sequences using novel feature development mechanisms. The obtained features were used to train the ensemble models, including blending, boosting, and bagging. The trained ensemble models were then subjected to independent testing and k-fold cross validation. The proposed model outperformed the preexisting predictors and revealed optimized scores based on major accuracy metrics."
    )
    #image = Image.open('pseudo.PNG')
    #st.image(image, width=400)

elif selected2 == "Predictor":
    #st.subheader("Predictor Page")
    import predictor

    exec(open('predictor.py').read())


elif selected2 == "Dataset":
    #st.subheader("Data Set")
    st.info("Positive Samples (modified m1A-sites)")
    with open("Sup1.fas", "rb") as file:
        btn = st.download_button(
            label="Download file",
            data=file,
            file_name="Sup1.fas",
            mime=""
        )
    st.info("Negative Samples (Non-modified m1A sites")
    with open("Sup2.fas", "rb") as file:
        btn = st.download_button(
            label="Download file",
            data=file,
            file_name="Sup2.fas",
            mime=""
        )