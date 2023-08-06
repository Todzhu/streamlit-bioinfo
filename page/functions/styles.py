import streamlit as st


@st.cache_data
def img_to_html(img_path):

    def img_to_bytes(img_path):
        import base64
        from pathlib import Path
        img_bytes = Path(img_path).read_bytes()
        encoded = base64.b64encode(img_bytes).decode()
        return encoded

    img_html = "<img src='data:image/png;base64,{}' style='width:100%;'>".format(img_to_bytes(img_path))
    
    return img_html

def add_color_to_cards():
    """
    Adds color to the expanders.
    Users don't need to call this function, is executed by default.
    """
    # Define your javascript
    my_js = """
                var cards = window.parent.document.getElementsByClassName("css-vhjbnf");
                for (var i = 0; i < cards.length; i++) {
                    let card = cards[i];
                    // See if there's content in the card
                    N_chars_in_cards = String(card.firstChild.innerHTML).length;
                    if (N_chars_in_cards >100){
                        card.style.border = "solid";
                        card.style.borderColor = "#20B2AA";
                        card.style.borderWidth = "2px";
                        card.style.padding = "10px";
                        card.style.borderRadius = "10px";
                        card.style.borderRadius = "10px";
                        card.addEventListener("mouseover", function(event){card.style.borderColor = "red"})
                        card.addEventListener("mouseout",  function(event){card.style.borderColor = "#20B2AA"})
                    }
                }    
            """
    # Wrapt the javascript as html code
    my_html = f"""
                <script>
                {my_js}
                </script>
                """
    # Execute your app
    st.components.v1.html(my_html, height=0, width=0)

    return
