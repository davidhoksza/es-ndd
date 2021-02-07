const path = require('path');
const webpack = require('webpack');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const { CleanWebpackPlugin } = require('clean-webpack-plugin');
const CopyPlugin = require("copy-webpack-plugin");


module.exports = {
  mode: 'development',
  entry: './src/js/index.js',
  devtool: 'source-map',
  devServer: {
    contentBase: './dist',
  },
  output: {
    filename: 'bundle.js',
    path: path.resolve(__dirname, 'dist'),
  },
  plugins: [
      new CleanWebpackPlugin(),
      new HtmlWebpackPlugin({
        template: 'src/index.html',
      }),
    new HtmlWebpackPlugin({
      template: 'src/contact.html',
      filename: "contact.html"
    }),
    new HtmlWebpackPlugin({
      template: './src/documentation.html',
      filename: "documentation.html"
    }),
    new webpack.DefinePlugin({
        'process.env.NODE_DEBUG': JSON.stringify(process.env.NODE_DEBUG)
      }),
      new CopyPlugin({
        patterns: [
          { from: path.resolve(__dirname, 'data/pdb'), to: path.resolve(__dirname, 'dist/data/pdb') },
          { from: path.resolve(__dirname, 'data/pdb-multi'), to: path.resolve(__dirname, 'dist/data/pdb-multi') },
          { from: path.resolve(__dirname, 'data/raptor'), to: path.resolve(__dirname, 'dist/data/raptor') },
          { from: path.resolve(__dirname, 'data/swissmodel'), to: path.resolve(__dirname, 'dist/data/swissmodel') },
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
    //   include: [path.resolve(__dirname, 'src')],
    //   loader: 'html-loader',
    // },
      {
      test: /.css$/,
      use: [{
        loader: "style-loader"
      }, {
        loader: "css-loader",
        options: {
          sourceMap: true
        }
      }]
    }, {
      test: /\.(woff(2)?|ttf|eot|jpg|svg)(\?v=\d+\.\d+\.\d+)?$/,
      type: 'asset/resource',
    }]
  },

}