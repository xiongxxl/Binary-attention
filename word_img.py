
import pandas as pd
from wordcloud import WordCloud
import matplotlib.pyplot as plt


def produce_word_image(dataframe):
    word_freq = dict(zip(dataframe['frag_smiles'], dataframe['Frequency']))
    wordcloud = WordCloud(width=800, height=400, background_color='white').generate_from_frequencies(word_freq)
    # 显示词云
    plt.figure(figsize=(10, 5))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis('off')
    plt.show()
    wordcloud.to_file('wordcloud_from_frequencies.png')

if __name__ == "__main__":

    data = {
        'frag_smiles': ['Python', 'Data', 'Science', 'Machine', 'Learning', 'AI', 'Python', 'Data', 'AI'],
        'Frequency': [10, 8, 6, 5, 4, 3, 10, 8, 3]
           }
    dataframe = pd.DataFrame(data)
    produce_word_image(dataframe)
