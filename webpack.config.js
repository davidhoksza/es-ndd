const path = require('path');
const webpack = require('webpack');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const { CleanWebpackPlugin } = require('clean-webpack-plugin');
const CopyPlugin = require("copy-webpack-plugin");

const TARGET_DIR = process.env.TARGET_DIR ? process.env.TARGET_DIR : 'DIST';

module.exports = {
  mode: process.env.MODE,
  entry: './src/js/index.js',
  devtool: 'source-map',
  devServer: {
    contentBase: `./${TARGET_DIR}`,
  },
  output: {
    filename: 'bundle.js',
    path: path.resolve(__dirname, TARGET_DIR),
  },
  plugins: [
      new CleanWebpackPlugin(),
      new HtmlWebpackPlugin({
        template: 'src/index.html',
      }),
    new HtmlWebpackPlugin({
      template: 'src/contact.html',
      filename: "contact.html",
      inject: false
    }),
    new HtmlWebpackPlugin({
      template: './src/documentation.html',
      filename: "documentation.html",
      inject: false
    }),
    new webpack.DefinePlugin({
        'process.env.NODE_DEBUG': JSON.stringify(process.env.NODE_DEBUG)
      }),
      new CopyPlugin({
        patterns: [
          { from: path.resolve(__dirname, 'data/v2/pdb'), to: path.resolve(__dirname, `${TARGET_DIR}/data/pdb`) },
          { from: path.resolve(__dirname, 'data/v2/pdb-multi'), to: path.resolve(__dirname, `${TARGET_DIR}/data/pdb-multi`) },
          { from: path.resolve(__dirname, 'data/v2/alphafold'), to: path.resolve(__dirname, `${TARGET_DIR}/data/alphafold`) },
          
          { from: path.resolve(__dirname, 'src/img'), to: path.resolve(__dirname, `${TARGET_DIR}/img`) },
        ],
      }),
  ],
  module: {
    rules: [{
      test: /\.(js|jsx)$/,
      include: [path.resolve(__dirname, 'src/js')],
      loader: 'babel-loader'
    },
    //   {
    //   test: /\.html$/i,
    //     use: {
    //       loader: 'html-loader',
    //     }
    //
    // },
    //   {
    //   test: /.css$/,
    //   use: [{
    //     loader: "style-loader"
    //   }, {
    //     loader: "css-loader",
    //     options: {
    //       sourceMap: true
    //     }
    //   }]
    // },
      {
        test: /\.css$/i,
        loader: 'file-loader',
        options: {
          outputPath: "css"
        }
      },
      {
      test: /\.(woff(2)?|ttf|eot|jpg|png|svg)(\?v=\d+\.\d+\.\d+)?$/,
      type: 'asset/resource',
    }]
  },

}